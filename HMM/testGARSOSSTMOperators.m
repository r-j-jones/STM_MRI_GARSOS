function results = testGARSOSSTMOperators(nufftRoot)
%testGARSOSSTMOPERATORS Run synthetic dimension and adjointness tests.
%
% Example:
%   results = testGARSOSSTMOperators( ...
%       fullfile(pwd, 'multi-scale-low-rank-MR-recon-master', 'nufft_toolbox'));

    if nargin > 0 && isfolder(nufftRoot)
        addpath(nufftRoot);
    end
    required = {'nufft_init', 'nufft', 'nufft_adj'};
    for i = 1:numel(required)
        if exist(required{i}, 'file') ~= 2
            error('testGARSOSSTMOperators:MissingDependency', ...
                'Add the NUFFT toolbox to the path before running tests.');
        end
    end

    nReadout = 8;
    nPartition = 4;
    nSpokes = 5;
    nbSpokesPerFrame = 2;
    nCoil = 3;
    nFrame = floor(nSpokes / nbSpokesPerFrame);
    nX = 12;
    nY = 10;
    nBasis = 3;

    sampled = complex(randn(nReadout*nPartition*nSpokes, nCoil), ...
        randn(nReadout*nPartition*nSpokes, nCoil));
    spokeValues = reshape(repmat(1:nSpokes, nReadout*nPartition, 1), [], 1);
    kx = repmat(linspace(-0.45, 0.45, nReadout).', 1, nSpokes);
    ky = sin((1:nSpokes) * pi / 5);
    ky = repmat(ky, nReadout, 1) .* cos(linspace(-0.45, 0.45, nReadout).');
    kxy = repmat(kx + 1i*ky, 1, 1, nFrame);
    coords4 = zeros(nReadout, nPartition, nSpokes, 3);
    coords4(:, :, :, 1) = permute(repmat(kx, 1, 1, nPartition), [1 3 2]);
    coords4(:, :, :, 2) = permute(repmat(ky, 1, 1, nPartition), [1 3 2]);
    kzValues = reshape(linspace(-0.375, 0.375, nPartition), [1 nPartition 1]);
    coords4(:, :, :, 3) = repmat(kzValues, nReadout, 1, nSpokes);
    coords = reshape(coords4, [], 3);

    prepared = prepareGARSOSSliceDataForMCNUFFT( ...
        sampled, coords, spokeValues, nReadout, nPartition, nSpokes, ...
        nbSpokesPerFrame, 2, ...
        'CoordinateTransform', diag([1 -1 -1]), ...
        'TranslationVector', [0; -1; -1]);
    expectedCoordinates = reshape(coords4, [], 3) * diag([1 -1 -1]);
    actualCoordinates = prepared.kSpaceLocations4D;
    actualCoordinates = reshape(actualCoordinates, [], 3);
    coordinateError = norm(actualCoordinates - expectedCoordinates) / ...
        max([norm(expectedCoordinates), eps]);
    if coordinateError > 1e-12
        error('testGARSOSSTMOperators:CoordinateTransform', ...
            'Coordinate transform did not preserve flattened sample order.');
    end
    validateGARSOSTrajectoryData(prepared.kxy, prepared.data);
    allZ = zeros(nReadout, nbSpokesPerFrame, nCoil, nFrame, nPartition, ...
        'like', sampled);
    for iz = 1:nPartition
        zPrepared = prepareGARSOSSliceDataForMCNUFFT( ...
            sampled, coords, spokeValues, nReadout, nPartition, nSpokes, ...
            nbSpokesPerFrame, iz);
        allZ(:, :, :, :, iz) = zPrepared.data;
    end
    data4 = reshape(sampled, [nReadout, nPartition, nSpokes, nCoil]);
    completeData4 = data4(:, :, 1:nFrame*nbSpokesPerFrame, :);
    kzRelativeError = abs(norm(allZ(:)) - norm(completeData4(:))) / ...
        max([norm(allZ(:)), norm(completeData4(:)), eps]);
    if prepared.nDiscardedSpokes ~= 1 || size(prepared.data, 2) ~= nbSpokesPerFrame
        error('testGARSOSSTMOperators:FrameGrouping', ...
            'Frame grouping or incomplete trailing-spoke handling failed.');
    end
    plans = buildGARSOSRadialNUFFTPlans(prepared.kxy, [nX nY]);
    senseMaps = complex(randn(nX, nY, nCoil), randn(nX, nY, nCoil));
    radial = buildGARSOSRadialEncodingOperator( ...
        prepared.kxy, [nX nY], senseMaps);
    ST = complex(randn(nX, nY, nFrame, nBasis), ...
        randn(nX, nY, nFrame, nBasis));
    stm = buildGARSOSSTMOperators(ST, senseMaps, radial);

    x = complex(randn(nX, nY, nFrame), randn(nX, nY, nFrame));
    y = complex(randn(nReadout, nbSpokesPerFrame, nCoil, nFrame), ...
        randn(nReadout, nbSpokesPerFrame, nCoil, nFrame));
    c = complex(randn(nX, nY, nBasis), randn(nX, nY, nBasis));
    radialX = radial.forward(x);
    radialAdjointY = radial.adjoint(y);
    stmC = stm.forward(c);
    stmAdjointY = stm.adjoint(y);
    lhsRadial = sum(conj(radialX(:)) .* y(:));
    rhsRadial = sum(conj(x(:)) .* radialAdjointY(:));
    lhsSTM = sum(conj(stmC(:)) .* y(:));
    rhsSTM = sum(conj(c(:)) .* stmAdjointY(:));
    dynamicX = stm.synthesis(c);
    synthesisAdjointX = stm.synthesisAdjoint(x);
    lhsSynthesis = sum(conj(dynamicX(:)) .* x(:));
    rhsSynthesis = sum(conj(c(:)) .* synthesisAdjointX(:));

    results = struct();
    results.kzTransformRelativeEnergyError = kzRelativeError;
    results.coordinateTransformRelativeError = coordinateError;
    results.radialAdjointRelativeError = abs(lhsRadial-rhsRadial) / ...
        max([abs(lhsRadial), abs(rhsRadial), eps]);
    results.stmSynthesisAdjointRelativeError = abs(lhsSynthesis-rhsSynthesis) / ...
        max([abs(lhsSynthesis), abs(rhsSynthesis), eps]);
    results.stmAdjointRelativeError = abs(lhsSTM-rhsSTM) / ...
        max([abs(lhsSTM), abs(rhsSTM), eps]);
    results.plans = plans;
    results.nFrame = nFrame;
    results.nDiscardedSpokes = prepared.nDiscardedSpokes;
    results.tolerance = 1e-8;
    if results.kzTransformRelativeEnergyError > results.tolerance || ...
            results.coordinateTransformRelativeError > results.tolerance || ...
            results.radialAdjointRelativeError > results.tolerance || ...
            results.stmSynthesisAdjointRelativeError > results.tolerance || ...
            results.stmAdjointRelativeError > results.tolerance
        error('testGARSOSSTMOperators:AdjointnessFailure', ...
            'One or more synthetic transform/adjointness tests exceeded tolerance.');
    end
    fprintf('kz transform relative energy error: %.3e\n', ...
        results.kzTransformRelativeEnergyError);
    fprintf('Radial adjoint relative error: %.3e\n', results.radialAdjointRelativeError);
    fprintf('STM synthesis adjoint relative error: %.3e\n', ...
        results.stmSynthesisAdjointRelativeError);
    fprintf('STM/radial adjoint relative error: %.3e\n', results.stmAdjointRelativeError);
end
