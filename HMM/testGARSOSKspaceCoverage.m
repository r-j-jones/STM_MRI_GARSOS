function results = testGARSOSKspaceCoverage()
%testGARSOSKSPACE COVERAGE Validate coverage metrics on synthetic trajectories.

    nSide = 16;
    [kx, ky] = ndgrid(linspace(-.45, .45, nSide));
    coordinates = [kx(:), ky(:), zeros(numel(kx), 1)];
    spokes = (1:numel(kx)).';
    dcf = ones(numel(kx), 1);
    uniform = analyzeGARSOSKspaceCoverage( ...
        coordinates, spokes, dcf, 1, numel(spokes), [16 16], [16 16]);
    if uniform.raw.CV > 1e-12 || uniform.raw.f_empty ~= 0
        error('testGARSOSKspaceCoverage:UniformCartesian', ...
            'Uniform Cartesian coverage did not produce uniform bins.');
    end

    hole = coordinates;
    hole(abs(hole(:, 1)) < .15 & abs(hole(:, 2)) < .15, :) = [];
    holeSpokes = (1:size(hole, 1)).';
    nonuniform = analyzeGARSOSKspaceCoverage( ...
        hole, holeSpokes, ones(size(hole, 1), 1), 1, size(hole, 1), ...
        [16 16], [16 16]);
    if nonuniform.raw.CV <= uniform.raw.CV || nonuniform.raw.f_empty <= 0
        error('testGARSOSKspaceCoverage:Hole', ...
            'A synthetic hole was not detected.');
    end

    nSpokes = 32;
    nReadout = 32;
    angles = (0:nSpokes - 1) * pi / nSpokes;
    radius = linspace(-.45, .45, nReadout);
    [radiusGrid, angleGrid] = ndgrid(radius, angles);
    radial = [radiusGrid(:) .* cos(angleGrid(:)), ...
        radiusGrid(:) .* sin(angleGrid(:)), zeros(numel(radiusGrid), 1)];
    radialSpokes = repmat((1:nSpokes), nReadout, 1);
    radialSpokes = radialSpokes(:);
    radialResult = analyzeGARSOSKspaceCoverage( ...
        radial, radialSpokes, ones(size(radial, 1), 1), nReadout, ...
        nSpokes, [16 16], [16 16]);
    if isempty(radialResult.angular) || ...
            any(~isfinite([radialResult.angular.CV]))
        error('testGARSOSKspaceCoverage:RadialAngular', ...
            'Radial angular coverage metrics are invalid.');
    end

    results = struct('uniformCV', uniform.raw.CV, ...
        'holeCV', nonuniform.raw.CV, ...
        'radialAngularCV', [radialResult.angular.CV]);
end
