function result = runGARSOSSTMIntegration(prepared, varargin)
%runGARSOSSTMINTEGRATION End-to-end GAR-SOS radial/STM operator assembly.
%
% Existing ST_maps, eigenValues, and kCal may be supplied; 
% otherwise kCal is generated and STM_computation is called when it is on the MATLAB path.

    p = inputParser;
    addParameter(p, 'ST_maps', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'eigenValues', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'kCal', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'nbSpokesPerFrame', [], ...
        @(x) isnumeric(x) && isscalar(x) && x > 0 && x == floor(x));
    addParameter(p, 'zSlice', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'CalibrationSize', [32 32], @(x) isnumeric(x) && numel(x)==2);
    addParameter(p, 'L', 7, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'STMOptions', {}, @(x) iscell(x));
    addParameter(p, 'RadialOptions', {}, @(x) iscell(x));
    addParameter(p, 'ReconstructionDir', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'ReferenceVolumeA', [], @(x) isempty(x) || isstruct(x));
    addParameter(p, 'ForceUpdate', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    if ~isstruct(prepared)
        error('runGARSOSSTMIntegration:InvalidPreparedData', 'prepared data or a manifest is required.');
    end
    if ~isfield(prepared, 'sampledData')
        if isempty(p.Results.ReferenceVolumeA) || isempty(p.Results.ReconstructionDir)
            error('runGARSOSSTMIntegration:MissingPreparationInputs', ...
                'A manifest requires ReconstructionDir and ReferenceVolumeA.');
        end
        prepared = prepareGARSOSSTMData(prepared, p.Results.ReconstructionDir, ...
            p.Results.ReferenceVolumeA, p.Results.ForceUpdate);
    end
    d = prepared.dimensions;
    if isempty(p.Results.nbSpokesPerFrame)
        error('runGARSOSSTMIntegration:MissingFrameSize', ...
            'nbSpokesPerFrame must be supplied explicitly.');
    end
    if isempty(p.Results.zSlice), zSlice = floor(d.imageSize(3)/2)+1; else, zSlice = p.Results.zSlice; end
    if ~isscalar(zSlice) || zSlice < 1 || zSlice ~= floor(zSlice) || ...
            zSlice > d.imageSize(3)
        error('runGARSOSSTMIntegration:InvalidSlice', ...
            'The direct slice-wise STM operator requires one valid zSlice.');
    end

    if isempty(p.Results.kCal)
        [kCal, prepOutputs, prepDiagnostics] = reconstructGARFrameForSTMCoordinateAware( ...
            prepared, p.Results.nbSpokesPerFrame, zSlice, p.Results.CalibrationSize);
    else
        kCal = p.Results.kCal; prepOutputs = struct(); prepDiagnostics = struct();
    end
    if isempty(p.Results.ST_maps)
        if isempty(which('STM_computation'))
            error('runGARSOSSTMIntegration:MissingSTM', ...
                'STM_computation is not on the MATLAB path; supply ST_maps or add STM_MRI_GARSOS-main.');
        end
        [ST_maps, eigenValues] = STM_computation(kCal, d.imageSize(1:2), ...
            p.Results.L, p.Results.STMOptions{:});
    else
        ST_maps = p.Results.ST_maps; eigenValues = p.Results.eigenValues;
    end
    if ndims(ST_maps) ~= 4 || size(ST_maps, 1) ~= d.imageSize(1) || ...
            size(ST_maps, 2) ~= d.imageSize(2)
        error('runGARSOSSTMIntegration:InvalidSTMMaps', ...
            'ST_maps must have dimensions [Nx Ny Nt L] matching the prepared grid.');
    end

    nFrame = size(ST_maps, 3);
    nReadout = prepared.nSamples;
    if nFrame ~= floor(prepared.nSpokes / p.Results.nbSpokesPerFrame)
        error('runGARSOSSTMIntegration:FrameMismatch', ...
            'ST_maps temporal dimension must match complete radial frames.');
    end
    preparedSlice = prepareGARSOSSliceDataForMCNUFFT( ...
        prepared.sampledData, prepared.kSpaceLocations, prepared.spokeIndex, ...
        prepared.nSamples, prepared.nPartitions, prepared.nSpokes, ...
        p.Results.nbSpokesPerFrame, zSlice, ...
        'CoordinateTransform', prepared.CA, ...
        'TranslationVector', prepared.Cb);
    validateGARSOSTrajectoryData(preparedSlice.kxy, preparedSlice.data);
    kxy = preparedSlice.kxy;
    radialData = preparedSlice.data;
    if size(radialData, 4) ~= nFrame
        error('runGARSOSSTMIntegration:RadialFrameMismatch', ...
            'Prepared radial data frame count does not match ST_maps.');
    end
    senseMapsSlice = reshape(prepared.sensitivityMaps(:, :, zSlice, :), ...
        [d.imageSize(1), d.imageSize(2), size(prepared.sensitivityMaps, 4)]);
    radialOperator = buildGARSOSRadialEncodingOperator(kxy, d.imageSize(1:2), ...
        senseMapsSlice, p.Results.RadialOptions{:});
    stmOperator = buildGARSOSSTMOperators(ST_maps, radialOperator.senseMaps, radialOperator);
    result = struct('kCal', kCal, 'ST_maps', ST_maps, 'eigenValues', eigenValues, ...
        'radialOperator', radialOperator, 'stmOperator', stmOperator, ...
        'radialData', radialData, 'trajectory', kxy, ...
        'senseMapsSlice', senseMapsSlice, ...
        'preparationOutputs', prepOutputs, 'diagnostics', prepDiagnostics, ...
        'nbSpokesPerFrame', p.Results.nbSpokesPerFrame);
end
