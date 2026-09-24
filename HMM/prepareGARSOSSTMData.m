function prepared = prepareGARSOSSTMData(metaManifest, reconstructionDir, referenceVolumeA, forceUpdate, varargin)
%prepareGARSOSSTMDATA Prepare shared GAR-SOS inputs for STM reconstruction.
%
% This is the common data boundary for GAR-SOS/STM experiments.  It uses
% the same radial preparation routine as radial_vibe_18_5_7_7_imFIAT_permuted
% and derives the coordinate transform from the reference geometry (rather
% than embedding an acquisition-specific transform).

    if nargin < 4 || isempty(forceUpdate), forceUpdate = false; end
    p = inputParser;
    addParameter(p, 'SequenceIndex', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
    addParameter(p, 'ImageSize', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
    addParameter(p, 'GridSize', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
    addParameter(p, 'KernelSize', [7 7 7], @(x) isnumeric(x) && numel(x) == 3);
    parse(p, varargin{:});

    if ~isstruct(metaManifest) || ~isfield(metaManifest, 'sequences')
        error('prepareGARSOSSTMData:InvalidManifest', ...
            'metaManifest.sequences is required.');
    end
    iSequence = p.Results.SequenceIndex;
    if iSequence > numel(metaManifest.sequences)
        error('prepareGARSOSSTMData:InvalidSequence', 'SequenceIndex is out of range.');
    end
    if ~isstruct(referenceVolumeA) || ~all(isfield(referenceVolumeA, {'R','v','r0','A'}))
        error('prepareGARSOSSTMData:InvalidReference', ...
            'referenceVolumeA must contain A, R, v, and r0.');
    end

    % Run the established radial pipeline so cached manifests and its input
    % conventions remain authoritative.
    [radialManifest, manifestFile] = radial_vibe_18_5_7_7_imFIAT_permuted( ...
        metaManifest, reconstructionDir, referenceVolumeA, forceUpdate, false);
    sequence = radialManifest.sequences(iSequence);
    if ~isfield(sequence, 'rawDataFile')
        sequence.rawDataFile = metaManifest.sequences(iSequence).rawDataFile;
    end

    s = readVD11MultiRaidFileStructure(sequence.rawDataFile, 'CalibrationOnly', true);
    q = mapVBVD(sequence.rawDataFile);
    if iscell(q), q = q{end}; end
    [noiseWhiteningTransform, channelSensitivityMaps, channelIds] = ...
        estimateCoilSensitivitieMaps2(s);
    [sampledData, kSpaceLocations, DCF, spokeIndex, timeStamps, nPartitions, ...
        nSpokes, nSamples] = prepareRadialVibeDataForGridding2_3( ...
        q, s, noiseWhiteningTransform, channelIds);

    referenceSize = size(referenceVolumeA.A);
    if isempty(p.Results.ImageSize)
        if referenceSize(3) == 64
            imageSize = [224 224 floor(referenceSize(3) * 1.5)];
        else
            imageSize = [224 224 96];
        end
    else
        imageSize = p.Results.ImageSize(:).';
    end
    if isempty(p.Results.GridSize)
        gridSize = ceil(1.375 * imageSize);
    else
        gridSize = p.Results.GridSize(:).';
    end

    dicomImageSize = referenceSize(1:3);
    visibleObjectSize = dicomImageSize;
    L1 = [eye(3), floor(dicomImageSize(:) / 2); 0 0 0 1];
    L2 = [diag([1 -1 -1]), [0; dicomImageSize(2)-1; dicomImageSize(3)-1]; 0 0 0 1];
    L3 = [referenceVolumeA.R * diag(referenceVolumeA.v), referenceVolumeA.r0(:); 0 0 0 1];
    Cs = L1 \ (L2 \ (L3 \ (eye(4) * L3 * L1)));
    CA = Cs(1:3,1:3);
    Cb = Cs(1:3,4);

    sensitivityMapArraySize = size(channelSensitivityMaps.A);
    periodicMaps = centeredArrayVolumeResize(channelSensitivityMaps, ...
        sensitivityMapArraySize(1:3) .* [3 1 3]);
    periodicMaps.A = repmat(channelSensitivityMaps.A, [3 1 3 1]);
    r0Machine = [s.m{2}.param.Meas{1}.DICOM.lGlobalTablePosSag, ...
        s.m{2}.param.Meas{1}.DICOM.lGlobalTablePosCor, ...
        s.m{2}.param.Meas{1}.DICOM.lGlobalTablePosTra];
    periodicMaps.r0 = periodicMaps.r0 + r0Machine(:);
    referenceResized = centeredArrayVolumeResize(referenceVolumeA, imageSize);
    alignedMaps = resampleArrayVolumeToSame(periodicMaps, referenceResized, 'linear', 0);

    prepared = struct('sampledData', sampledData, 'kSpaceLocations', kSpaceLocations, ...
        'DCF', DCF, 'densityCompensation', DCF, 'spokeIndex', spokeIndex, ...
        'iiSpokes', spokeIndex, 'timeStamps', timeStamps, ...
        'nPartitions', nPartitions, 'nSpokes', nSpokes, 'nSamples', nSamples, ...
        'dimensions', struct('imageSize', imageSize, 'gridSize', gridSize, ...
            'dicomImageSize', dicomImageSize, 'visibleObjectSize', visibleObjectSize, ...
            'nCoils', size(sampledData,2)), ...
        'sensitivityMaps', alignedMaps.A, 'alignedSensitivityMaps', alignedMaps.A, ...
        'channelSensitivityMaps', alignedMaps, 'Cs', Cs, 'CA', CA, 'Cb', Cb, ...
        'L1', L1, 'L2', L2, 'L3', L3, 'manifest', radialManifest, ...
        'manifestFile', manifestFile, 'sequenceIndex', iSequence, ...
        'kernelSize', p.Results.KernelSize);
end
