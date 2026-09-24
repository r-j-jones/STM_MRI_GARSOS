function op = buildGARSOSSTMOperators(ST_maps, senseMaps, radialOperator)
%buildGARSOSSTMOPERATORS Compose STM synthesis with radial encoding.
%
% ST_maps:  [Nx, Ny, Nt, L]
% senseMaps: [Nx, Ny, Nc]
% C:         [Nx, Ny, L]
% dynamic:   [Nx, Ny, Nt]

    validateattributes(ST_maps, {'numeric'}, {'4d', 'nonempty'});
    validateattributes(senseMaps, {'numeric'}, {'3d', 'nonempty'});
    requiredFields = {'imageSize', 'nFrame', 'nCoil', 'forward', 'adjoint'};
    if ~all(isfield(radialOperator, requiredFields))
        error('buildGARSOSSTMOperators:InvalidRadialOperator', ...
            'radialOperator is missing required fields.');
    end

    [nX, nY, nFrame, nBasis] = size(ST_maps);
    if size(senseMaps, 1) ~= nX || size(senseMaps, 2) ~= nY || ...
            size(senseMaps, 3) ~= radialOperator.nCoil || ...
            nFrame ~= radialOperator.nFrame
        error('buildGARSOSSTMOperators:DimensionMismatch', ...
            'STM, sensitivity-map, and radial-operator dimensions disagree.');
    end

    op = struct();
    op.imageSize = [nX nY];
    op.nFrame = nFrame;
    op.nBasis = nBasis;
    op.nCoil = size(senseMaps, 3);
    op.ST_maps = ST_maps;
    op.senseMaps = senseMaps;
    op.radialOperator = radialOperator;
    op.synthesis = @synthesis;
    op.synthesisAdjoint = @synthesisAdjoint;
    op.forward = @(C) radialOperator.forward(synthesis(C));
    op.adjoint = @(y) synthesisAdjoint(radialOperator.adjoint(y));
    op.forwardVector = @(c) reshape(op.forward(reshape(c, [nX nY nBasis])), [], 1);
    op.adjointVector = @(y) reshape(op.adjoint( ...
        reshape(y, [radialOperator.nReadout, radialOperator.nSpokes, ...
        radialOperator.nCoil, nFrame])), [], 1);
    op.normalVector = @(c) op.adjointVector(op.forwardVector(c));

    function X = synthesis(C)
        validateattributes(C, {'numeric'}, {'size', [nX nY nBasis]});
        X = sum(ST_maps .* reshape(C, [nX nY 1 nBasis]), 4);
    end

    function C = synthesisAdjoint(X)
        validateattributes(X, {'numeric'}, {'size', [nX nY nFrame]});
        C = sum(conj(ST_maps) .* reshape(X, [nX nY nFrame 1]), 3);
    end
end
