function result = main_deid_permuted_STM_GARSOS(inputPath, outputPath, opts)
%main_deid_permuted_STM_GARSOS Prepare and reconstruct GAR-SOS STM data.
%
% This is a new GAR-SOS-specific entry point based on the preparation
% sequence in main_deid_permuted_STM. Existing production entry points and
% STM source files are intentionally left unchanged.
%
% Required:
%   inputPath, outputPath, opts.nbSpokesPerFrame
%
% Optional opts fields:
%   forceUpdate, referenceVolumeA, ST_maps, eigenValues, kCal,
%   zSlice, calibrationSize, L, STMOptions, RadialOptions

    % % Add path to Fessler toolbox for NUFFT operations
    % addpath('/path/to/multi-scale-low-rank-MR-recon-master/nufft_toolbox');

    narginchk(3, 3);
    if ~isstruct(opts) || ~isfield(opts, 'nbSpokesPerFrame')
        error('main_deid_permuted_STM_GARSOS:MissingFrameSize', ...
            'opts.nbSpokesPerFrame is required.');
    end
    validateattributes(opts.nbSpokesPerFrame, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    if ~isfolder(inputPath) || ~isfolder(outputPath)
        error('main_deid_permuted_STM_GARSOS:InvalidPath', ...
            'inputPath and outputPath must be existing directories.');
    end

    if ~isfield(opts, 'forceUpdate'), opts.forceUpdate = false; end
    if ~isfield(opts, 'zSlice'), opts.zSlice = []; end
    if ~isfield(opts, 'calibrationSize'), opts.calibrationSize = [32 32]; end
    if ~isfield(opts, 'L'), opts.L = 7; end
    if ~isfield(opts, 'STMOptions'), opts.STMOptions = {}; end
    if ~isfield(opts, 'RadialOptions'), opts.RadialOptions = {}; end

    reconstructionDir = fullfile(inputPath, 'reconstruction');
    inputFile = fullfile(reconstructionDir, 'input.mat');
    if exist(inputFile, 'file') ~= 2
        error('main_deid_permuted_STM_GARSOS:MissingInput', ...
            'Expected reconstruction input file: %s', inputFile);
    end
    loaded = load(inputFile, 'kSpaceID1', 'kSpaceID2');
    if ~isfield(loaded, 'kSpaceID1') || ~isfield(loaded, 'kSpaceID2')
        error('main_deid_permuted_STM_GARSOS:MissingKSpaceIds', ...
            'input.mat must contain kSpaceID1 and kSpaceID2.');
    end

    if ~isfield(opts, 'referenceVolumeA') || isempty(opts.referenceVolumeA)
        [referenceVolumeA, ~, ~] = loadReferenceDicomImage(inputPath, true);
    else
        referenceVolumeA = opts.referenceVolumeA;
    end
    [kSpaceID, ~] = parseKspaceIds(loaded.kSpaceID1, loaded.kSpaceID2);
    [metaManifest, ~] = radial_vibe_18_5_4_timestamp_imFIAT( ...
        kSpaceID, inputPath, outputPath, opts.forceUpdate);

    prepared = prepareGARSOSSTMData( ...
        metaManifest, outputPath, referenceVolumeA, opts.forceUpdate);

    integrationArgs = { ...
        'nbSpokesPerFrame', opts.nbSpokesPerFrame, ...
        'zSlice', opts.zSlice, ...
        'CalibrationSize', opts.calibrationSize, ...
        'L', opts.L, ...
        'STMOptions', opts.STMOptions, ...
        'RadialOptions', opts.RadialOptions};
    for name = {'ST_maps', 'eigenValues', 'kCal'}
        if isfield(opts, name{1}) && ~isempty(opts.(name{1}))
            integrationArgs(end + 1:end + 2) = {name{1}, opts.(name{1})}; %#ok<AGROW>
        end
    end
    result = runGARSOSSTMIntegration(prepared, integrationArgs{:});
    result.prepared = prepared;
    result.inputPath = inputPath;
    result.outputPath = outputPath;
end
