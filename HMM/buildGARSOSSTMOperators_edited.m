function op = buildGARSOSSTMOperators(ST_maps, senseMaps, radialOperator)
%buildGARSOSSTMOperators Compose STM synthesis with radial encoding.
%
% ST_maps:   [Nx, Ny, Nt, L]
% senseMaps: [Nx, Ny, Nc]
% C:         [Nx, Ny, L]
% dynamic:   [Nx, Ny, Nt]

    validateattributes(ST_maps, {'numeric'}, {'nonempty'});
    validateattributes(senseMaps, {'numeric'}, {'nonempty'});

    requiredFields = { ...
        'imageSize', 'nFrame', 'nCoil', ...
        'nReadout', 'nSpokes', 'forward', 'adjoint'};

    if ~all(isfield(radialOperator, requiredFields))
        error('buildGARSOSSTMOperators:InvalidRadialOperator', ...
            'radialOperator is missing required fields.');
    end

    [nX, nY, nFrame, nBasis] = size(ST_maps);

    if ~isequal(radialOperator.imageSize, [nX nY]) || ...
            size(senseMaps, 1) ~= nX || ...
            size(senseMaps, 2) ~= nY || ...
            size(senseMaps, 3) ~= radialOperator.nCoil || ...
            nFrame ~= radialOperator.nFrame

        error('buildGARSOSSTMOperators:DimensionMismatch', ...
            ['STM maps, sensitivity maps, and radial-operator ' ...
             'dimensions disagree.']);
    end

    nCoil = size(senseMaps, 3);
    nReadout = radialOperator.nReadout;
    nSpokes = radialOperator.nSpokes;

    % Cache these handles locally.
    radialForward = radialOperator.forward;
    radialAdjoint = radialOperator.adjoint;

    % Avoid recomputing this in every adjoint operation.
    conjugateSTMaps = conj(ST_maps);

    op = struct();
    op.imageSize = [nX nY];
    op.nFrame = nFrame;
    op.nBasis = nBasis;
    op.nCoil = nCoil;
    op.ST_maps = ST_maps;
    op.senseMaps = senseMaps;
    op.radialOperator = radialOperator;

    op.synthesis = @synthesis;
    op.synthesisAdjoint = @synthesisAdjoint;
    op.forward = @forward;
    op.adjoint = @adjoint;
    op.forwardVector = @forwardVector;
    op.adjointVector = @adjointVector;
    op.normalVector = @normalVector;

    function X = synthesis(C)
        % C: [nX, nY, nBasis]
        % X: [nX, nY, nFrame]
        X = sum(ST_maps .* reshape(C, [nX nY 1 nBasis]), 4);
    end

    function C = synthesisAdjoint(X)
        % X: [nX, nY, nFrame]
        % C: [nX, nY, nBasis]
        C = sum(conjugateSTMaps .* ...
            reshape(X, [nX nY nFrame 1]), 3);
    end

    function y = forward(C)
        X = synthesis(C);
        y = radialForward(X);
    end

    function C = adjoint(y)
        X = radialAdjoint(y);
        C = synthesisAdjoint(X);
    end

    function yVector = forwardVector(cVector)
        C = reshape(cVector, [nX nY nBasis]);
        y = forward(C);
        yVector = y(:);
    end

    function cVector = adjointVector(yVector)
        y = reshape(yVector, ...
            [nReadout, nSpokes, nCoil, nFrame]);

        C = adjoint(y);
        cVector = C(:);
    end

    function cVector = normalVector(cVector)
        cVector = adjointVector(forwardVector(cVector));
    end
end