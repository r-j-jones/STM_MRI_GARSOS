function summary = sweepGARSOSKspaceCoverage(prepared, spokesPerFrameValues, ...
    calibrationSizes, varargin)
%sweepGARSOSKSPACECOVERAGE Compare coverage across frame and square sizes.

    p = inputParser;
    addParameter(p, 'FrameIndex', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NBins', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
    parse(p, varargin{:});

    if ~isstruct(prepared) || ~isfield(prepared, 'kSpaceLocations') || ...
            ~isfield(prepared, 'spokeIndex') || ~isfield(prepared, 'DCF')
        error('sweepGARSOSKspaceCoverage:InvalidPreparedData', ...
            'prepared must contain kSpaceLocations, spokeIndex, and DCF.');
    end
    spokesPerFrameValues = spokesPerFrameValues(:);
    calibrationSizes = reshape(calibrationSizes, [], 2);
    imageSize = prepared.dimensions.imageSize(1:2);
    rows = cell(numel(spokesPerFrameValues) * size(calibrationSizes, 1), 1);
    rowIndex = 0;
    for iFrameSize = 1:numel(spokesPerFrameValues)
        for iSquare = 1:size(calibrationSizes, 1)
            rowIndex = rowIndex + 1;
            squareSize = calibrationSizes(iSquare, :);
            options = {'FrameIndex', p.Results.FrameIndex};
            if ~isempty(p.Results.NBins)
                options = [options, {'NBins', p.Results.NBins}]; %#ok<AGROW>
            end
            metrics = analyzeGARSOSKspaceCoverage( ...
                prepared.kSpaceLocations, prepared.spokeIndex, prepared.DCF, ...
                prepared.nSamples, spokesPerFrameValues(iFrameSize), ...
                squareSize, imageSize, options{:});
            rows{rowIndex} = {spokesPerFrameValues(iFrameSize), squareSize(1), ...
                squareSize(2), metrics.raw.CV, metrics.dcf.CV, ...
                metrics.raw.f_empty, metrics.raw.f_low25, ...
                metrics.dcf.f_empty, metrics.dcf.f_low25, ...
                metrics.nearestNeighbor.CV, ...
                localAngularCV(metrics.angular)};
        end
    end
    summary = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'spokesPerFrame', 'calibrationNx', 'calibrationNy', ...
        'rawCV', 'dcfCV', 'rawEmptyFraction', 'rawLow25Fraction', ...
        'dcfEmptyFraction', 'dcfLow25Fraction', ...
        'nearestNeighborCV', 'angularGapCV'});
end

function value = localAngularCV(metrics)
    values = [metrics.CV];
    values = values(isfinite(values));
    if isempty(values)
        value = NaN;
    else
        value = mean(values);
    end
end
