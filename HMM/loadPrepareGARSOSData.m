function prepared = loadPrepareGARSOSData( ...
    metaManifest, referenceVolumeA, griddingParams, varargin)
%loadPrepareGARSOSDATA Load and prepare one GAR-SOS sequence for reuse.
%
% The returned struct contains the raw/preprocessed radial arrays, aligned
% sensitivity maps, image/grid dimensions, and the geometry transforms used
% by both calibration and direct STM reconstruction.

    p = inputParser;
    addParameter(p, 'SequenceIndex', 1, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
    addParameter(p, 'MakePlots', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    iSequence = p.Results.SequenceIndex;
    makePlots = p.Results.MakePlots;

    if ~isstruct(metaManifest) || ~isfield(metaManifest, 'sequences')
        error('loadPrepareGARSOSData:InvalidManifest', ...
            'metaManifest.sequences is required.');
    end
    if iSequence > numel(metaManifest.sequences)
        error('loadPrepareGARSOSData:InvalidSequence', ...
            'SequenceIndex is out of range.');
    end
    if numel(metaManifest.sequences) ~= 1
        warning('loadPrepareGARSOSData:SingleSequenceRequired', ...
            'GAR-SOS STM preparation supports exactly one sequence.');
    end

    requiredNames = {'R', 'v', 'r0', 'A'};
    validStruct = isstruct(referenceVolumeA) && ...
        all(cellfun(@(name) isfield(referenceVolumeA, name), requiredNames));    
    validArrayVolume = isa(referenceVolumeA, 'ArrayVolume') && ...
        all(cellfun(@(name) isprop(referenceVolumeA, name), requiredNames));    
    if ~(validStruct || validArrayVolume)
        error('loadPrepareGARSOSData:InvalidReference', ...
            ['referenceVolumeA must be either a struct or an ArrayVolume ' ...
             'containing R, v, r0, and A.']);
    end

    sequence = metaManifest.sequences(iSequence);
    rawDataFile = sequence.rawDataFile;
    fprintf('\n-- Preparing GAR-SOS sequence %d/%d --\n', ...
        iSequence, numel(metaManifest.sequences));

    s = readVD11MultiRaidFileStructure(rawDataFile, 'CalibrationOnly', true);
    q = mapVBVD(rawDataFile);
    if iscell(q)
        q = q{end};
    end

    [noiseWhiteningTransform, channelSensitivityMaps, channelIds, debugStructure] = ...
        estimateGARSOSCoilSensitivityMaps(s);
    % [orig_noiseWhiteningTransform, orig_channelSensitivityMaps, orig_channelIds, ~] = ...
    %     estimateCoilSensitivitieMaps2(s);

    
    if makePlots
        plotDir = fullfile(metaManifest.dstRootDir,'coil-plots');
        if ~exist(plotDir,'dir'), mkdir(plotDir); end
    
        % Plot noise covariance + whitening matrices
        figsNoise = plotNoiseMatrices(debugStructure);
        fnames = fieldnames(figsNoise);
        for ff=1:length(fnames)
            curname = fnames{ff};
            curfig = figsNoise.(curname);
            print(curfig,fullfile(plotDir,sprintf('%s.png',curname)),'-dpng','-r300');
        end
    
        % Plot coil ref images
        figsCoils = plotCoilSensitivityDiagnostics(debugStructure, [65,72,130]);
        fnames = fieldnames(figsCoils);
        for ff=1:length(fnames)
            curname = fnames{ff};
            curfig = figsCoils.(curname);
            print(curfig,fullfile(plotDir,sprintf('%s.png',curname)),'-dpng','-r300');
        end
    end



    [sampledData, kSpaceLocations, densityCompensation, iiSpokes, ...
        timeStamps, nPartitions, nSpokes, nSamples] = ...
        prepareRadialVibeDataForGridding2_3( ...
        q, s, noiseWhiteningTransform, channelIds);

    referenceSize = size(referenceVolumeA.A);
    if ~isfield(griddingParams, 'imageSize') || isempty(griddingParams.imageSize)
        if referenceSize(3) == 64
            imageSize = [224 224 floor(referenceSize(3) * 1.5)];
        else
            imageSize = [224 224 96];
        end
    else
        imageSize = griddingParams.imageSize(:).';
    end
    if ~isfield(griddingParams, 'gridOversamplingFactor') || ...
            isempty(griddingParams.gridOversamplingFactor)
        griddingParams.gridOversamplingFactor = 1.375;
    end
    if ~isfield(griddingParams, 'VisibleObjectSize') || ...
            isempty(griddingParams.VisibleObjectSize)
        griddingParams.VisibleObjectSize = referenceSize(1:3);
    end
    gridSize = ceil(griddingParams.gridOversamplingFactor .* imageSize);

    dicomImageSize = referenceSize(1:3);
    visibleObjectSize = dicomImageSize;
    L1 = [eye(3), floor(dicomImageSize(:) / 2); 0 0 0 1];
    L2 = [diag([1 -1 -1]), ...
        [0; dicomImageSize(2)-1; dicomImageSize(3)-1]; 0 0 0 1];
    L3 = [referenceVolumeA.R * diag(referenceVolumeA.v), ...
        referenceVolumeA.r0(:); 0 0 0 1];
    Cs = L1 \ (L2 \ (L3 \ (eye(4) * L3 * L1)));
    CA = Cs(1:3, 1:3);
    Cb = Cs(1:3, 4);

    sensitivityMapArraySize = size(channelSensitivityMaps.A);
    periodicMaps = centeredArrayVolumeResize( ...
        channelSensitivityMaps, sensitivityMapArraySize(1:3) .* [3 1 3]);
    periodicMaps.A = repmat(channelSensitivityMaps.A, [3 1 3 1]);
    r0Machine = [ ...
        s.m{2}.param.Meas{1}.DICOM.lGlobalTablePosSag, ...
        s.m{2}.param.Meas{1}.DICOM.lGlobalTablePosCor, ...
        s.m{2}.param.Meas{1}.DICOM.lGlobalTablePosTra];
    periodicMaps.r0 = periodicMaps.r0 + r0Machine(:);
    referenceVolumeC = centeredArrayVolumeResize(referenceVolumeA, imageSize);
    alignedMaps = resampleArrayVolumeToSame( ...
        periodicMaps, referenceVolumeC, 'linear', 0);

    if makePlots
        % Plot aligned coil sens maps (1 slice, all coils)
        coilmapPlotSlices = [25,50,75];
        for pp=1:length(coilmapPlotSlices)
            coilmapPlotSlice = coilmapPlotSlices(pp);
    
            figSensMaps = figure;
            im('row',4,'col',5,abs(squeeze(alignedMaps.A(:,:,coilmapPlotSlice,:))));
            title(sprintf('Coil sensitivity maps (z=%d)',coilmapPlotSlice),'FontSize',20);
            cbb=colorbar;
            cbb.FontSize = 15;
            figSensMaps.Position = [710 127 1065 765];
            figSensMaps.Color = 'w';
            drawnow;
    
            outfigname = fullfile(plotDir, sprintf('CoilSensMaps_z%d.png',coilmapPlotSlice));
            print(gcf,outfigname,'-dpng','-r300');
    
        end
    end

    prepared = struct();
    prepared.sequenceIndex = iSequence;
    prepared.rawDataFile = rawDataFile;
    prepared.sampledData = sampledData;
    prepared.kSpaceLocations = kSpaceLocations;
    prepared.densityCompensation = densityCompensation;
    prepared.DCF = densityCompensation;
    prepared.spokeIndex = iiSpokes;
    prepared.iiSpokes = iiSpokes;
    prepared.timeStamps = timeStamps;
    prepared.nPartitions = nPartitions;
    prepared.nSpokes = nSpokes;
    prepared.nSamples = nSamples;
    prepared.nChannels = size(sampledData, 2);
    prepared.sensitivityMaps = alignedMaps.A;
    prepared.alignedSensitivityMaps = alignedMaps.A;
    prepared.channelSensitivityMaps = alignedMaps;
    prepared.Cs = Cs;
    prepared.CA = CA;
    prepared.Cb = Cb;
    prepared.L1 = L1;
    prepared.L2 = L2;
    prepared.L3 = L3;
    prepared.dimensions = struct( ...
        'imageSize', imageSize, ...
        'gridSize', gridSize, ...
        'referenceSize', referenceSize, ...
        'dicomImageSize', dicomImageSize, ...
        'visibleObjectSize', visibleObjectSize, ...
        'nCoils', size(sampledData, 2));
    prepared.metadata = struct( ...
        'sensitivityMapArraySize', sensitivityMapArraySize, ...
        'r0Machine', r0Machine, ...
        'referenceSize', referenceSize);
    prepared.griddingParams = griddingParams;
    prepared.raw = struct( ...
        'sampledData', sampledData, ...
        'kSpaceLocations', kSpaceLocations, ...
        'densityCompensation', densityCompensation, ...
        'spokeIndex', iiSpokes, ...
        'timeStamps', timeStamps, ...
        'nPartitions', nPartitions, ...
        'nSpokes', nSpokes, ...
        'nSamples', nSamples);
    prepared.geometry = struct( ...
        'Cs', Cs, 'CA', CA, 'Cb', Cb, ...
        'L1', L1, 'L2', L2, 'L3', L3,...
        'visibleObjectSize', visibleObjectSize);
    prepared.coil = struct( ...
        'sensitivityMaps', alignedMaps.A, ...
        'alignedSensitivityMaps', alignedMaps.A, ...
        'channelSensitivityMaps', alignedMaps, ...
        'nChannels', size(sampledData, 2), ...
        'arraySize', sensitivityMapArraySize, ...
        'r0Machine', r0Machine);
end
