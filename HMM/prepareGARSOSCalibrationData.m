function [manifest, manifestFile, stmCalibOutputs] = ...
    prepareGARSOSCalibrationData( ...
    metaManifest, reconstructionDir, referenceVolumeA, frameParams, ...
    griddingParams, forceUpdate, varargin)
%prepareGARSOSCALIBRATIONDATA Prepare and save GAR-SOS STM calibration data.
%
% Required inputs:
%   metaManifest       GAR-SOS sequence manifest.
%   reconstructionDir  Destination reconstruction directory.
%   referenceVolumeA   Reference ArrayVolume used for geometry alignment.
%   frameParams        Struct with nSpokesPerFrame, zSlice, calibrationSize.
%   griddingParams     Struct controlling image/grid reconstruction. Missing
%                      fields receive the defaults used by the prior code.
%
% Optional inputs:
%   forceUpdate        Defaults to false. May be logical or a numeric
%                      integer equal to 0 or 1. The value
%                      'RerunWithoutSaving' forces recomputation without
%                      saving the resulting manifest.
%
% Outputs:
%   manifest           Saved manifest metadata.
%   manifestFile       Path to the saved manifest.mat file.
%   stmCalibOutputs    Large reconstruction arrays for the current run:
%                      frameData, reconstructedVolumes, sensitivitySlice,
%                      senseImages, and kSpace2D. These are intentionally
%                      not saved in the manifest. If a cached manifest is
%                      used, this output is returned as an empty struct.
%
% Large stmCalibOutputs data are unavailable after a cached-manifest load.
% 
% Rerun with forceUpdate=true to regenerate and return them. Use
%  forceUpdate='RerunWithoutSaving' to regenerate without saving.
% 
% Robert Jones  |  09-24-2026


    if nargin == 0
        localDisplayUsage();
        manifest = [];
        manifestFile = [];
        stmCalibOutputs = [];
        return;
    end
    
    narginchk(4, inf);

    if nargin < 5 || isempty(griddingParams)
        griddingParams = struct();
    end
    if nargin < 6 || isempty(forceUpdate)
        forceUpdate = false;
    end
    preparedData = [];
    if ~isempty(varargin)
        if mod(numel(varargin), 2) ~= 0
            error('prepareGARSOSCalibrationData:InvalidOptions', ...
                'Optional arguments must be name/value pairs.');
        end
        for iOption = 1:2:numel(varargin)
            if strcmpi(varargin{iOption}, 'PreparedData')
                preparedData = varargin{iOption + 1};
            else
                error('prepareGARSOSCalibrationData:UnknownOption', ...
                    'Unknown optional argument.');
            end
        end
    end

    rerunWithoutSaving = ...
        (ischar(forceUpdate) && isrow(forceUpdate) && ...
            strcmp(forceUpdate, 'RerunWithoutSaving')) || ...
        (isstring(forceUpdate) && isscalar(forceUpdate) && ...
            forceUpdate == "RerunWithoutSaving");

    if rerunWithoutSaving
        forceUpdate = true;
    elseif islogical(forceUpdate) || isnumeric(forceUpdate)
        validateattributes(forceUpdate, {'logical', 'numeric'}, ...
            {'scalar', 'real', 'finite', 'integer', '>=', 0, '<=', 1});
        forceUpdate = logical(forceUpdate);
    else
        error('prepareGARSOSCalibrationData:InvalidForceUpdate', ...
            ['forceUpdate must be logical, a numeric integer equal to 0 or 1, ' ...
             'or ''RerunWithoutSaving''.']);
    end

    frameParams = localValidateFrameParams(frameParams);
    validateattributes(metaManifest, {'struct'}, {'scalar'});
    validateattributes(referenceVolumeA, {'struct'}, {'scalar'});
    if ~isfield(metaManifest, 'sequences') || ...
            ~isstruct(metaManifest.sequences) || ...
            numel(metaManifest.sequences) ~= 1
        error('prepareGARSOSCalibrationData:SingleSequenceRequired', ...
            'GAR-SOS STM preparation supports exactly one sequence.');
    end

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
        iSequence = 1;
        if ~isfield(manifest, 'sequences') || ...
                ~isstruct(manifest.sequences) || ...
                numel(manifest.sequences) ~= 1 || ...
                ~isfield(manifest.sequences(1), 'calibData') || ...
                ~isstruct(manifest.sequences(1).calibData) || ...
                ~isscalar(manifest.sequences(1).calibData)
            error('prepareGARSOSCalibrationData:MissingCachedCalibData', ...
                'Cached manifest must contain one sequence with scalar calibData.');
        end
        cached = manifest.sequences(1).calibData;
        if ~isfield(cached, 'kCal') || ~isnumeric(cached.kCal) || ...
                isempty(cached.kCal) || ~isfield(cached, 'nbSpokesPerFrame') || ...
                ~isnumeric(cached.nbSpokesPerFrame) || ...
                ~isscalar(cached.nbSpokesPerFrame) || ...
                cached.nbSpokesPerFrame < 1 || ...
                cached.nbSpokesPerFrame ~= floor(cached.nbSpokesPerFrame)
            error('prepareGARSOSCalibrationData:InvalidCachedCalibData', ...
                'Cached calibData.kCal and calibData.nbSpokesPerFrame are invalid.');
        end
        if ~isfield(cached, 'zSlice') || ~isnumeric(cached.zSlice) || ...
                isempty(cached.zSlice) || any(cached.zSlice < 1) || ...
                any(cached.zSlice ~= floor(cached.zSlice)) || ...
                ~isfield(cached, 'calibrationSize') || ...
                ~isnumeric(cached.calibrationSize) || ...
                numel(cached.calibrationSize) ~= 2 || ...
                any(cached.calibrationSize < 1) || ...
                any(cached.calibrationSize ~= floor(cached.calibrationSize))
            error('prepareGARSOSCalibrationData:InvalidCachedCalibData', ...
                'Cached calibData.zSlice and calibData.calibrationSize are invalid.');
        end
        return;
    end

    permissiveMakeDirIfNotExist(destinationDir);
    nSequences = 1;
    manifest = metaManifest;
    totalTimer = tic;

    for iSequence = 1:nSequences
        sequenceTimer = tic;
        sequence = metaManifest.sequences(iSequence);
        if isempty(preparedData)
            preparedData = loadPrepareGARSOSData( ...
                metaManifest, referenceVolumeA, griddingParams, ...
                'SequenceIndex', iSequence);
        end
        rawDataFile = preparedData.rawDataFile;
        sampledData = preparedData.sampledData;
        kSpaceLocations = preparedData.kSpaceLocations;
        densityCompensation = preparedData.densityCompensation;
        iiSpokes = preparedData.spokeIndex;
        timeStamps = preparedData.timeStamps;
        nPartitions = preparedData.nPartitions;
        nSpokes = preparedData.nSpokes;
        nSamples = preparedData.nSamples;

        % Set dimensions
        referenceSize = preparedData.dimensions.referenceSize;
        imageSize = preparedData.dimensions.imageSize;
        gridSize = preparedData.dimensions.gridSize;
        nChannels = preparedData.nChannels;
        dicomImageSize = preparedData.dimensions.dicomImageSize;
        visibleObjectSize = preparedData.dimensions.visibleObjectSize;

        % Derive the existing physical-coordinate transform; do not
        % hard-code acquisition-specific CA, Cb, or Cs values.
        L1 = preparedData.L1;
        L2 = preparedData.L2;
        L3 = preparedData.L3;
        Cs = preparedData.Cs;
        CA = preparedData.CA;
        Cb = preparedData.Cb;

        % Prepare aligned coil sensitivity maps
        sensitivityMapArraySize = preparedData.metadata.sensitivityMapArraySize;
        r0Machine = preparedData.metadata.r0Machine;
        alignedMaps = preparedData.channelSensitivityMaps;

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
        manifest.sequences(iSequence).kCal = kCal;
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
    if ~isempty(frameParams.zSlice)
        validateattributes(frameParams.zSlice, {'numeric'}, ...
            {'vector', 'integer', 'positive'});
    end
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

% prepareGARSOSCalibrationData( ...
%     metaManifest, reconstructionDir, referenceVolumeA, frameParams, ...
%     griddingParams, forceUpdate)

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

    fprintf('%% Force update, griddingParams omitted:\n');
    fprintf([ ...
        'prepareGARSOSCalibrationData(' ...
        'meta, dir, ref, frameParams, [], true)\n\n']);

    fprintf('%% Force update and explicit griddingParams:\n');
    fprintf([ ...
        'prepareGARSOSCalibrationData(' ...
        'meta, dir, ref, frameParams, griddingParams, true)\n\n']);

    fprintf('%% Recompute without saving:\n');
    fprintf([ ...
        'prepareGARSOSCalibrationData( ...\n' ...
        '    meta, dir, ref, ' ...
        'frameParams, griddingParams, ''RerunWithoutSaving'')\n\n']);

    fprintf('forceUpdate values:\n');
    fprintf('  false or 0              Use a current cached manifest.\n');
    fprintf('  true or 1               Recompute and save the manifest.\n');
    fprintf([ ...
        '  ''RerunWithoutSaving''  Recompute without saving ' ...
        'the manifest.\n\n']);
end