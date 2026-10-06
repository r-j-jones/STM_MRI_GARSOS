function results = analyzeGARSOSKspaceCoverage( ...
    sampleCoordinates, spokeIndex, densityCompensation, ...
    nReadout, nSpokesPerFrame, calibrationSize, imageSize, varargin)
%analyzeGARSOSKSPACE COVERAGE Quantify trajectory coverage in a calibration square.
%
% Coordinates are the normalized cycles/FOV coordinates used by
% prepareRadialVibeDataForGridding2_3. The calibration square follows the
% same centered-index convention as reconstructGARFrameForSTM: its half
% width is calibrationSize ./ (2 .* imageSize(1:2)).

    p = inputParser;
    addParameter(p, 'FrameIndex', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
    addParameter(p, 'NBins', calibrationSize, @(x) isnumeric(x) && numel(x) == 2 && all(x >= 1 & x == floor(x)));
    addParameter(p, 'RadialBands', [0 .25 .5 .75 1], @(x) isnumeric(x) && isvector(x) && numel(x) >= 2 && x(1) == 0 && x(end) == 1 && all(diff(x) > 0));
    addParameter(p, 'DCFZeroTolerance', 1e-12, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'UseTwoPiAngles', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    validateattributes(sampleCoordinates, {'numeric'}, {'2d', 'ncols', 3});
    validateattributes(spokeIndex, {'numeric'}, {'vector', 'numel', size(sampleCoordinates, 1)});
    validateattributes(densityCompensation, {'numeric'}, {'vector', 'numel', size(sampleCoordinates, 1)});
    validateattributes(nReadout, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(nSpokesPerFrame, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(calibrationSize, {'numeric'}, {'vector', 'numel', 2, 'integer', 'positive'});
    validateattributes(imageSize, {'numeric'}, {'vector', 'numel', 2, 'positive'});

    nRows = size(sampleCoordinates, 1);
    if mod(nRows, nReadout) ~= 0
        error('analyzeGARSOSKspaceCoverage:InvalidReadoutGrouping', ...
            'The sample count must be divisible by nReadout.');
    end
    if any(~isfinite(sampleCoordinates(:))) || any(~isfinite(spokeIndex(:))) || ...
            any(~isfinite(densityCompensation(:)))
        error('analyzeGARSOSKspaceCoverage:NonFiniteInput', ...
            'Coordinates, spoke identifiers, and DCF values must be finite.');
    end

    spokeValues = unique(spokeIndex(:), 'stable');
    nFrames = floor(numel(spokeValues) / nSpokesPerFrame);
    if p.Results.FrameIndex > nFrames
        error('analyzeGARSOSKspace:InvalidFrame', ...
            'FrameIndex exceeds the number of complete frames.');
    end
    frameSpokes = spokeValues((p.Results.FrameIndex - 1) * nSpokesPerFrame + ...
        (1:nSpokesPerFrame));
    inFrame = ismember(spokeIndex, frameSpokes);

    kxy = sampleCoordinates(:, 1:2);
    halfWidth = calibrationSize(:).' ./ (2 .* imageSize(:).');
    square = all(abs(kxy) <= halfWidth, 2) & inFrame;
    valid = square & all(isfinite(kxy), 2);
    kxy = kxy(valid, :);
    dcf = densityCompensation(valid);
    spokes = spokeIndex(valid);

    nBins = p.Results.NBins(:).';
    edgesX = linspace(-halfWidth(1), halfWidth(1), nBins(1) + 1);
    edgesY = linspace(-halfWidth(2), halfWidth(2), nBins(2) + 1);
    [ix, iy, binValid] = localBinIndices(kxy, edgesX, edgesY);
    dcf = dcf(binValid);
    kxy = kxy(binValid, :);
    spokes = spokes(binValid);
    H = accumarray([iy ix], 1, [nBins(2) nBins(1)], @sum, 0);
    Hw = accumarray([iy ix], dcf, [nBins(2) nBins(1)], @sum, 0);
    Hw(abs(Hw) <= p.Results.DCFZeroTolerance) = 0;
    validBins = true(size(H));
    results = struct();
    results.frameIndex = p.Results.FrameIndex;
    results.frameSpokes = frameSpokes;
    results.coordinateSystem = 'normalized cycles/FOV';
    results.calibrationSize = calibrationSize;
    results.calibrationHalfWidth = halfWidth;
    results.binEdges = {edgesX, edgesY};
    results.raw = localDensityMetrics(H, validBins, 'raw');
    results.dcf = localDensityMetrics(Hw, validBins, 'dcf');
    results.raw.map = H;
    results.raw.normalizedMap = localNormalize(H, validBins);
    results.dcf.map = Hw;
    results.dcf.normalizedMap = localNormalize(Hw, validBins);
    results.sampleCoordinates = kxy;
    results.sampleWeights = dcf;

    results.nearestNeighbor = localNearestNeighbor(kxy);
    results.angular = localAngularMetrics(kxy, spokes, halfWidth, ...
        p.Results.RadialBands, p.Results.UseTwoPiAngles);
end

function [ix, iy, valid] = localBinIndices(points, edgesX, edgesY)
    ix = discretize(points(:, 1), edgesX);
    iy = discretize(points(:, 2), edgesY);
    valid = isfinite(ix) & isfinite(iy);
    ix = ix(valid);
    iy = iy(valid);
end

function metrics = localDensityMetrics(H, validBins, name)
    values = H(validBins);
    meanValue = mean(values);
    metrics = struct('name', name, 'mean', meanValue, ...
        'std', std(values), 'CV', localCV(values), ...
        'CV_percent', 100 * localCV(values), ...
        'f_empty', mean(values == 0), ...
        'f_empty_percent', 100 * mean(values == 0), ...
        'f_low25', mean(values < .25 * meanValue), ...
        'f_low25_percent', 100 * mean(values < .25 * meanValue), ...
        'f_low50', mean(values < .50 * meanValue), ...
        'f_low50_percent', 100 * mean(values < .50 * meanValue));
end

function normalized = localNormalize(H, validBins)
    meanValue = mean(H(validBins));
    if meanValue > 0
        normalized = H ./ meanValue;
    else
        normalized = zeros(size(H), 'like', H);
    end
end

function value = localCV(values)
    meanValue = mean(values);
    if meanValue == 0
        value = NaN;
    else
        value = std(values) / meanValue;
    end
end

function metrics = localNearestNeighbor(points)
    nPoints = size(points, 1);
    metrics = struct('distances', [], 'mean', NaN, 'std', NaN, ...
        'CV', NaN, 'percentiles', [NaN NaN NaN]);
    if nPoints < 2
        return;
    end
    if exist('knnsearch', 'file') == 2
        [~, distances] = knnsearch(points, points, 'K', 2);
        distances = distances(:, 2);
    else
        % Toolbox-free fallback; this is O(N^2) in the selected square.
        distances = inf(nPoints, 1);
        for iPoint = 1:nPoints
            delta = points - points(iPoint, :);
            delta(iPoint, :) = [];
            distances(iPoint) = min(sqrt(sum(delta.^2, 2)));
        end
    end
    metrics.distances = distances;
    metrics.mean = mean(distances);
    metrics.std = std(distances);
    metrics.CV = metrics.std / max(metrics.mean, eps);
    metrics.percentiles = localPercentiles(distances, [5 50 95]);
end

function metrics = localAngularMetrics(kxy, spokeIndex, halfWidth, bands, useTwoPi)
    radius = sqrt(sum((kxy ./ halfWidth).^2, 2));
    period = 2*pi;
    if ~useTwoPi
        period = pi;
    end
    angle = mod(atan2(kxy(:, 2), kxy(:, 1)), period);
    metrics = repmat(struct('radialRange', [], 'angles', [], ...
        'gaps', [], 'meanGap', NaN, 'stdGap', NaN, 'CV', NaN, ...
        'minGap', NaN, 'maxGap', NaN), 1, numel(bands) - 1);
    for iBand = 1:numel(metrics)
        selected = radius >= bands(iBand) & radius <= bands(iBand + 1);
        selectedAngles = angle(selected);
        selectedSpokes = spokeIndex(selected);
        [~, firstIndex] = unique(selectedSpokes, 'stable');
        selectedAngles = sort(selectedAngles(firstIndex));
        if numel(selectedAngles) >= 2
            gaps = diff([selectedAngles; selectedAngles(1) + period]);
            meanGap = mean(gaps);
            stdGap = std(gaps);
        else
            gaps = [];
            meanGap = NaN;
            stdGap = NaN;
        end
        metrics(iBand).radialRange = bands(iBand:iBand + 1);
        metrics(iBand).angles = selectedAngles;
        metrics(iBand).gaps = gaps;
        metrics(iBand).meanGap = meanGap;
        metrics(iBand).stdGap = stdGap;
        metrics(iBand).CV = stdGap / max(meanGap, eps);
        if ~isempty(gaps)
            metrics(iBand).minGap = min(gaps);
            metrics(iBand).maxGap = max(gaps);
        end
    end
end

function values = localPercentiles(data, percentages)
    data = sort(data(isfinite(data)));
    if isempty(data)
        values = NaN(size(percentages));
        return;
    end
    positions = 1 + (numel(data) - 1) .* percentages ./ 100;
    lower = floor(positions);
    upper = ceil(positions);
    fraction = positions - lower;
    values = data(lower) + fraction .* (data(upper) - data(lower));
end
