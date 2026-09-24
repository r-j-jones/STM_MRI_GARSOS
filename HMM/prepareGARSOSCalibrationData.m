function [manifest, manifestFile, calibData, stmCalibOutputs] = ...
    prepareGARSOSCalibrationData( ...
    metaManifest, reconstructionDir, referenceVolumeA, forceUpdate, ...
    frameParams, griddingParams)
%prepareGARSOSCALIBRATIONDATA Prepare and save GAR-SOS STM calibration data.
%
% Required inputs:
%   metaManifest       GAR-SOS sequence manifest.
%   reconstructionDir  Destination reconstruction directory.
%   referenceVolumeA   Reference ArrayVolume used for geometry alignment.
%   frameParams        Struct with nSpokesPerFrame, zSlice, calibrationSize.
%
% Optional inputs:
%   forceUpdate        Defaults to false. May be logical or a numeric
%                      integer equal to 0 or 1. The value
%                      'ReruunWithoutSaving' forces recomputation without
%                      saving the resulting manifest.
%   griddingParams     Struct controlling image/grid reconstruction. Missing
%                      fields receive the defaults used by the prior code.
%
% Outputs:
%   manifest           Saved manifest metadata.
%   manifestFile       Path to the saved manifest.mat file.
%   calibData          Calibration data for the current run, including kCal, 
%                      nSpokesPerFrame, zSlice, and calibrationSize.
%   stmCalibOutputs    Large reconstruction arrays for the current run:
%                      frameData, reconstructedVolumes, sensitivitySlice,
%                      senseImages, and kSpace2D. These are intentionally
%                      not saved in the manifest. If a cached manifest is
%                      used, this output is returned as an empty struct.
%
% Large stmCalibOutputs data are unavailable after a cached-manifest load.
% 
% Rerun with forceUpdate=true to regenerate and return them. Use
%  forceUpdate='ReruunWithoutSaving' to regenerate without saving.
% 
% Robert Jones  |  09-24-2026


    if nargin == 0
        localDisplayUsage();
        manifest = [];
        manifestFile = [];
        calibData = [];
        stmCalibOutputs = [];
        return;
    end
    
    narginchk(4, 6);

    % Because forceUpdate precedes frameParams in the function signature,
    % shift the arguments when forceUpdate is omitted.
    if nargin <= 5 && isstruct(forceUpdate)
        if nargin == 5
            griddingParams = frameParams;
        else
            griddingParams = struct();
        end
        frameParams = forceUpdate;
        forceUpdate = false;
    elseif nargin < 6 || isempty(griddingParams)
        griddingParams = struct();
    end

    rerunWithoutSaving = ...
        (ischar(forceUpdate) && isrow(forceUpdate) && ...
            strcmp(forceUpdate, 'ReruunWithoutSaving')) || ...
        (isstring(forceUpdate) && isscalar(forceUpdate) && ...
            forceUpdate == "ReruunWithoutSaving");

    if rerunWithoutSaving
        forceUpdate = true;
    elseif islogical(forceUpdate) || isnumeric(forceUpdate)
        validateattributes(forceUpdate, {'logical', 'numeric'}, ...
            {'scalar', 'real', 'finite', 'integer', '>=', 0, '<=', 1});
        forceUpdate = logical(forceUpdate);
    else
        error('prepareGARSOSCalibrationData:InvalidForceUpdate', ...
            ['forceUpdate must be logical, a numeric integer equal to 0 or 1, ' ...
             'or ''ReruunWithoutSaving''.']);
    end

    frameParams = localValidateFrameParams(frameParams);
    validateattributes(metaManifest, {'struct'}, {'scalar'});
    validateattributes(referenceVolumeA, {'struct'}, {'scalar'});

    griddingParams = localNormalizeGriddingParams( ...
        griddingParams, referenceVolumeA, frameParams);
    stmCalibOutputs = struct();

    destinationDir = fullfile( ...
        reconstructionDir, sprintf('GARSOS_STM_Calib'));
    manifestFile = fullfile(destinationDir, 'garsos-stm-calib-manifest.mat');
    reconstructionTimestamp = '2018-04-18 13:58:25 -04:00';

    [manifestLoaded, manifest, ~] = ...
        loadIfExistAndTimestampNewerThan(manifestFile, reconstructionTimestamp);
    if manifestLoaded && ~forceUpdate
        fprintf('Using cached GAR-SOS calibration manifest: %s\n', manifestFile);
        return;
    end

    permissiveMakeDirIfNotExist(destinationDir);
    nSequences = numel(metaManifest.sequences);
    manifest = metaManifest;
    totalTimer = tic;

    for iSequence = 1:nSequences
        sequenceTimer = tic;
        sequence = metaManifest.sequences(iSequence);
        rawDataFile = sequence.rawDataFile;

        fprintf('\n-- Preparing GAR-SOS sequence %d/%d --\n', ...
            iSequence, nSequences);
        s = readVD11MultiRaidFileStructure( ...
            rawDataFile, 'CalibrationOnly', true);
        q = mapVBVD(rawDataFile);
        if iscell(q)
            q = q{end};
        end

        % Preprocessing
        [noiseWhiteningTransform, channelSensitivityMaps, channelIds, ~] = ...
            estimateCoilSensitivitieMaps2(s);
        [sampledData, kSpaceLocations, densityCompensation, iiSpokes, ...
            timeStamps, nPartitions, nSpokes, nSamples] = ...
            prepareRadialVibeDataForGridding2_3( ...
            q, s, noiseWhiteningTransform, channelIds);

        % Set dimensions
        referenceSize = size(referenceVolumeA.A);
        imageSize = griddingParams.imageSize;
        gridSize = ceil(griddingParams.gridOversamplingFactor .* imageSize);
        nChannels = size(sampledData, 2);
        dicomImageSize = referenceSize;
        visibleObjectSize = griddingParams.VisibleObjectSize;

        % Derive the existing physical-coordinate transform; do not
        % hard-code acquisition-specific CA, Cb, or Cs values.
        L1 = [eye(3), floor(dicomImageSize(:) / 2); 0 0 0 1];
        L2 = [diag([1 -1 -1]), ...
            [0; dicomImageSize(2) - 1; dicomImageSize(3) - 1]; 0 0 0 1];
        L3 = [referenceVolumeA.R * diag(referenceVolumeA.v), ...
            referenceVolumeA.r0(:); 0 0 0 1];
        Cs = L1 \ (L2 \ (L3 \ (eye(4) * L3 * L1)));
        CA = Cs(1:3, 1:3);
        Cb = Cs(1:3, 4);

        % Prepare aligned coil sensitivity maps
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

        % Store small sequence metadata that are not represented elsewhere.
        sequenceMetadata = struct();
        sequenceMetadata.iiSpokes = 0:nSpokes - 1;
        sequenceMetadata.complexFile = fullfile( ...
            destinationDir, sprintf('Sequence_%06.0f_ComplexData.raw', ...
            sequence.sequenceId));
        sequenceMetadata.imageSize = imageSize;
        sequenceMetadata.dicomImageSize = dicomImageSize;
        sequenceMetadata.nSpokes = nSpokes;
        sequenceMetadata.nPartitions = nPartitions;
        sequenceMetadata.nSamples = nSamples;
        sequenceMetadata.nChannels = nChannels;
        sequenceMetadata.timeStamps = timeStamps;
        sequenceMetadata.sensitivityMapArraySize = sensitivityMapArraySize;
        sequenceMetadata.Cs = Cs;
        sequenceMetadata.CA = CA;
        sequenceMetadata.Cb = Cb;

        % Structs with other params for generating STM calibration data
        currentFrameParams = frameParams;
        currentFrameParams.zSlice = localResolveZSlice( ...
            frameParams.zSlice, imageSize(3));
        currentFrameParams.calibrationSize = frameParams.calibrationSize;

        currentGriddingParams = griddingParams;
        currentGriddingParams.gridSize = gridSize;
        currentGriddingParams.nDimensions = 3;
        currentGriddingParams.referenceSize = referenceSize;
        currentGriddingParams.dicomImageSize = dicomImageSize;
        currentGriddingParams.nChannels = nChannels;
        currentGriddingParams.r0_machine = r0Machine;

        inputOptions = localBuildInputOptions( ...
            currentFrameParams, currentGriddingParams, imageSize, ...
            gridSize, CA, Cb, visibleObjectSize);

        % Generate STM calibration data
        fprintf('--  Generating STM calibration data..\n');
        fprintf(' Spokes per frame = %d,\n Calibration size = %d x %d\n', ...
            currentFrameParams.nSpokesPerFrame,...
            currentFrameParams.calibrationSize);

        [kCal, outputs, diagnostics] = reconstructGARFrameForSTM( ...
            sampledData, kSpaceLocations, densityCompensation, iiSpokes, ...
            alignedMaps, imageSize, gridSize, ...
            currentFrameParams.nSpokesPerFrame, ...
            currentFrameParams.zSlice, ...
            currentFrameParams.calibrationSize, ...
            'CoordinateTransform', inputOptions.CoordinateTransform, ...
            'TranslationVector', inputOptions.TranslationVector, ...
            'VisibleObjectSize', inputOptions.VisibleObjectSize, ...
            'PhaseWidth', inputOptions.PhaseWidth, ...
            'SpatialSigma', inputOptions.SpatialSigma, ...
            'DensityScale', inputOptions.DensityScale, ...
            'UseGpu', inputOptions.UseGpu, ...
            'DisplaySlice', inputOptions.DisplaySlice, ...
            'DisplaySliceIndex', inputOptions.DisplaySliceIndex, ...
            'ReturnDiagnostics', inputOptions.ReturnDiagnostics, ...
            'KernelSize', inputOptions.KernelSize, ...
            'OutputIsGpuArray', inputOptions.OutputIsGpuArray);

        % Keep large arrays available to the caller, but never save them.
        if isfield(diagnostics, 'frameData')
            stmCalibOutputs.frameData = diagnostics.frameData;
            diagnostics = rmfield(diagnostics, 'frameData');
        else
            stmCalibOutputs.frameData = [];
        end
        stmCalibOutputs.reconstructedVolumes = outputs.reconstructedVolumes;
        stmCalibOutputs.sensitivitySlice = outputs.sensitivitySlice;
        stmCalibOutputs.senseImages = outputs.senseImages;
        stmCalibOutputs.kSpace2D = outputs.kSpace2D;
        diagnostics.timers = outputs.timers;

        % Easy access to STM calibration data variables of interest
        calibData = struct( ...
            'kCal', kCal, ...
            'nbSpokesPerFrame', currentFrameParams.nSpokesPerFrame, ...
            'zSlice', currentFrameParams.zSlice, ...
            'calibrationSize', currentFrameParams.calibrationSize);

        % Store results in manifest struct
        manifest.sequences(iSequence) = localMergeStruct( ...
            manifest.sequences(iSequence), sequenceMetadata);
        manifest.sequences(iSequence).diagnostics = diagnostics;
        manifest.sequences(iSequence).elapProcessingTime = toc(sequenceTimer);
        manifest.sequences(iSequence).calibData = calibData;
        manifest.sequences(iSequence).inputOptions = inputOptions;
        manifest.sequences(iSequence).griddingParams = currentGriddingParams;
        manifest.sequences(iSequence).frameParams = currentFrameParams;
    end

    % Save manifest struct to mat file
     manifest.totalElapsedTime = toc(totalTimer);
    if ~rerunWithoutSaving
        saveWithTimestamp(manifestFile, manifest, iso8601Now());
        [manifestLoaded, manifest, ~] = ...
            loadIfExistAndTimestampNewerThan( ...
                manifestFile, reconstructionTimestamp);
        assert(manifestLoaded, ...
            'prepareGARSOSCalibrationData:ManifestSaveFailed', ...
            'The calibration manifest could not be reloaded after saving.');
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    LOCAL HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frameParams = localValidateFrameParams(frameParams)
%localValidateFrameParams Validate frameParams

    if ~isstruct(frameParams) || ~isscalar(frameParams)
        error('prepareGARSOSCalibrationData:InvalidFrameParams', ...
            'frameParams must be a scalar struct.');
    end
    required = {'nSpokesPerFrame', 'zSlice', 'calibrationSize'};
    if ~all(isfield(frameParams, required))
        error('prepareGARSOSCalibrationData:MissingFrameParams', ...
            'frameParams must contain nSpokesPerFrame, zSlice, calibrationSize.');
    end
    validateattributes(frameParams.nSpokesPerFrame, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(frameParams.zSlice, {'numeric'}, ...
        {'vector', 'integer', 'positive'});
    validateattributes(frameParams.calibrationSize, {'numeric'}, ...
        {'vector', 'numel', 2, 'integer', 'positive'});
end

function zSlice = localResolveZSlice(zSlice, nZ)
    if isempty(zSlice)
        zSlice = 1:nZ;
    elseif any(zSlice > nZ)
        error('prepareGARSOSCalibrationData:InvalidZSlice', ...
            'frameParams.zSlice exceeds the reconstructed image volume.');
    else
        zSlice = unique(zSlice(:).', 'stable');
    end
end

function params = localNormalizeGriddingParams(params, referenceVolumeA, frameParams)
%localNormalizeGriddingParams Set default griddingParams

    if ~isstruct(params) || ~isscalar(params)
        error('prepareGARSOSCalibrationData:InvalidGriddingParams', ...
            'griddingParams must be a scalar struct.');
    end
    referenceSize = size(referenceVolumeA.A);
    if referenceSize(3) == 64
        defaultImageSize = [224 224 floor(referenceSize(3) * 1.5)];
    else
        defaultImageSize = [224 224 96];
    end
    defaults = struct( ...
        'gridOversamplingFactor', 1.375, ...
        'kernelSize', [7 7 7], ...
        'spatialSigma', 0, ...
        'imageSize', defaultImageSize, ...
        'PhaseWidth', 1, ...
        'UseGpu', true, ...
        'DisplaySlice', true, ...
        'DisplaySliceIndex', [], ...
        'ReturnDiagnostics', true, ...
        'DensityScale', 1, ...
        'OutputIsGpuArray', false, ...
        'VisibleObjectSize', referenceSize);
    fields = fieldnames(defaults);
    for iField = 1:numel(fields)
        name = fields{iField};
        if ~isfield(params, name) || isempty(params.(name))
            params.(name) = defaults.(name);
        end
    end
    validateattributes(params.gridOversamplingFactor, {'numeric'}, ...
        {'scalar', 'positive', 'finite'});
    validateattributes(params.kernelSize, {'numeric'}, ...
        {'vector', 'numel', 3, 'positive', 'integer'});
    validateattributes(params.spatialSigma, {'numeric'}, ...
        {'scalar', 'nonnegative', 'finite'});
    validateattributes(params.imageSize, {'numeric'}, ...
        {'vector', 'numel', 3, 'positive', 'integer'});
    validateattributes(params.PhaseWidth, {'numeric'}, ...
        {'scalar', 'nonzero', 'finite'});
    validateattributes(params.VisibleObjectSize, {'numeric'}, ...
        {'vector', 'numel', 3, 'positive', 'finite'});
    validateattributes(params.UseGpu, {'logical'}, {'scalar'});
    validateattributes(params.DisplaySlice, {'logical'}, {'scalar'});
    validateattributes(params.ReturnDiagnostics, {'logical'}, {'scalar'});
    validateattributes(params.OutputIsGpuArray, {'logical'}, {'scalar'});
    if isempty(params.DisplaySliceIndex)
        params.DisplaySliceIndex = floor(params.imageSize(3) / 2) + 1;
    end
    if params.DisplaySliceIndex < 1 || ...
            params.DisplaySliceIndex > params.imageSize(3)
        error('prepareGARSOSCalibrationData:InvalidDisplaySlice', ...
            'DisplaySliceIndex must be within imageSize(3).');
    end
    if any(frameParams.calibrationSize > params.imageSize(1:2))
        error('prepareGARSOSCalibrationData:InvalidCalibrationSize', ...
            'calibrationSize cannot exceed imageSize in-plane dimensions.');
    end
end

function inputOptions = localBuildInputOptions( ...
    frameParams, griddingParams, imageSize, gridSize, CA, Cb, visibleObjectSize)
%localBuildInputOptions Store variables in structs

    inputOptions = struct();
    inputOptions.imageSize = imageSize;
    inputOptions.gridSize = gridSize;
    inputOptions.nSpokesPerFrame = frameParams.nSpokesPerFrame;
    inputOptions.zSlice = frameParams.zSlice;
    inputOptions.calibrationSize = frameParams.calibrationSize;
    inputOptions.CoordinateTransform = CA;
    inputOptions.TranslationVector = Cb;
    inputOptions.VisibleObjectSize = visibleObjectSize;
    inputOptions.PhaseWidth = griddingParams.PhaseWidth;
    inputOptions.SpatialSigma = griddingParams.spatialSigma;
    inputOptions.DensityScale = griddingParams.DensityScale;
    inputOptions.UseGpu = griddingParams.UseGpu;
    inputOptions.DisplaySlice = griddingParams.DisplaySlice;
    inputOptions.DisplaySliceIndex = griddingParams.DisplaySliceIndex;
    inputOptions.ReturnDiagnostics = griddingParams.ReturnDiagnostics;
    inputOptions.KernelSize = griddingParams.kernelSize;
    inputOptions.OutputIsGpuArray = griddingParams.OutputIsGpuArray;
end

function merged = localMergeStruct(original, additions)
%localMergeStruct Merge struct content

    merged = original;
    names = fieldnames(additions);
    for iName = 1:numel(names)
        merged.(names{iName}) = additions.(names{iName});
    end
end

function localDisplayUsage()
%localDisplayUsage Display documentation and supported call forms.

    fprintf('\n%s\n', help(mfilename));

    fprintf('Supported call forms:\n\n');

    fprintf('%% forceUpdate omitted, griddingParams omitted:\n');
    fprintf([ ...
        'prepareGARSOSCalibrationData(' ...
        'meta, dir, ref, frameParams)\n\n']);

    fprintf('%% forceUpdate omitted, griddingParams supplied:\n');
    fprintf([ ...
        'prepareGARSOSCalibrationData(' ...
        'meta, dir, ref, frameParams, griddingParams)\n\n']);

    fprintf('%% Existing argument order, griddingParams omitted:\n');
    fprintf([ ...
        'prepareGARSOSCalibrationData(' ...
        'meta, dir, ref, true, frameParams)\n\n']);

    fprintf('%% Existing full argument order:\n');
    fprintf([ ...
        'prepareGARSOSCalibrationData(' ...
        'meta, dir, ref, 1, frameParams, griddingParams)\n\n']);

    fprintf('%% Recompute without saving:\n');
    fprintf([ ...
        'prepareGARSOSCalibrationData( ...\n' ...
        '    meta, dir, ref, ''ReruunWithoutSaving'', ' ...
        'frameParams, griddingParams)\n\n']);

    fprintf('forceUpdate values:\n');
    fprintf('  false or 0              Use a current cached manifest.\n');
    fprintf('  true or 1               Recompute and save the manifest.\n');
    fprintf([ ...
        '  ''ReruunWithoutSaving''  Recompute without saving ' ...
        'the manifest.\n\n']);
end