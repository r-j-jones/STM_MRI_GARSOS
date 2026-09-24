function report = validateGARSOSTrajectoryData(kxy, data, varargin)
%validateGARSOSTrajectoryData Verify flattened trajectory/data correspondence.
%
% kxy  [nReadout, nSpokes, nFrame], normalized complex kx+i*ky
% data [nReadout, nSpokes, nCoil, nFrame]

    p = inputParser;
    addParameter(p, 'Tolerance', 1e-12, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, varargin{:});

    validateattributes(kxy, {'numeric'}, {'3d', 'nonempty'});
    validateattributes(data, {'numeric'}, {'4d', 'nonempty'});
    if size(kxy, 1) ~= size(data, 1) || size(kxy, 2) ~= size(data, 2) || ...
            size(kxy, 3) ~= size(data, 4)
        error('validateGARSOSTrajectoryData:SizeMismatch', ...
            'Trajectory and data dimensions do not share readout, spoke, and frame axes.');
    end

    trajectoryFlat = kxy(:);
    dataFlat = reshape(data, size(data, 1) * size(data, 2), size(data, 3), size(data, 4));
    dataFlat = dataFlat(:);
    if isempty(trajectoryFlat) || isempty(dataFlat)
        error('validateGARSOSTrajectoryData:EmptyInput', ...
            'Trajectory and data must be nonempty.');
    end

    report = struct();
    report.trajectoryShape = size(kxy);
    report.dataShape = size(data);
    report.flattenedTrajectoryLength = numel(trajectoryFlat);
    report.flattenedDataLength = numel(dataFlat);
    report.flatteningOrder = 'readout fastest, spoke next, frame outermost; coil is retained as the middle data axis';
    report.lengthsAgree = report.flattenedTrajectoryLength * size(data, 3) == numel(data);
    report.finiteTrajectory = all(isfinite(kxy(:)));
    report.finiteData = all(isfinite(data(:)));
    report.tolerance = p.Results.Tolerance;

    if ~report.lengthsAgree || ~report.finiteTrajectory || ~report.finiteData
        error('validateGARSOSTrajectoryData:InvalidInput', ...
            'Trajectory/data correspondence validation failed.');
    end
end
