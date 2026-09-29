function op = buildGARSOSSTMOperators( ...
    ST_maps, senseMaps, radialOperator, varargin)
%buildGARSOSSTMOperators Compose STM synthesis with radial encoding.
%
% ST_maps:   [Nx, Ny, Nt, L]
% senseMaps: [Nx, Ny, Nc]
% C:         [Nx, Ny, L]
% dynamic:   [Nx, Ny, Nt]
%
% Optional name-value argument:
%   'UseSinglePrecision'  Convert STM maps, sensitivity maps, inputs, and
%                         outputs to single precision. Default: false.
%
% If omitted and radialOperator.useSinglePrecision exists, the STM
% operator inherits that setting.
%
% Examples:
%   op = buildGARSOSSTMOperators( ...
%       ST_maps, senseMaps, radialOperator);
%
%   op = buildGARSOSSTMOperators( ...
%       ST_maps, senseMaps, radialOperator, ...
%       'UseSinglePrecision', true);

    validateattributes(ST_maps, {'numeric'}, {'nonempty'});
    validateattributes(senseMaps, {'numeric'}, {'nonempty'});

    p = inputParser;
    p.FunctionName = mfilename;

    if isfield(radialOperator, 'useSinglePrecision')
        defaultUseSingle = radialOperator.useSinglePrecision;
    else
        defaultUseSingle = false;
    end

    addParameter(p, 'UseSinglePrecision', defaultUseSingle, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x) && ...
        (x == 0 || x == 1));

    parse(p, varargin{:});
    useSinglePrecision = logical(p.Results.UseSinglePrecision);

    requiredFields = { ...
        'imageSize', 'nFrame', 'nCoil', ...
        'nReadout', 'nSpokes', 'forward', 'adjoint'};

    if ~all(isfield(radialOperator, requiredFields))
        error('buildGARSOSSTMOperators:InvalidRadialOperator', ...
            'radialOperator is missing required fields.');
    end

    if isfield(radialOperator, 'useSinglePrecision') && ...
            logical(radialOperator.useSinglePrecision) ~= ...
            useSinglePrecision
        error('buildGARSOSSTMOperators:PrecisionMismatch', ...
            ['STM and radial operators must use the same precision. ' ...
             'Rebuild both with matching UseSinglePrecision values.']);
    end

    if useSinglePrecision
        ST_maps = single(ST_maps);
        senseMaps = single(senseMaps);
        precisionClass = 'single';
    else
        precisionClass = class(ST_maps);
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

    radialForward = radialOperator.forward;
    radialAdjoint = radialOperator.adjoint;

    conjugateSTMaps = conj(ST_maps);

    op = struct();
    op.imageSize = [nX nY];
    op.nFrame = nFrame;
    op.nBasis = nBasis;
    op.nCoil = nCoil;
    op.ST_maps = ST_maps;
    op.senseMaps = senseMaps;
    op.radialOperator = radialOperator;
    op.useSinglePrecision = useSinglePrecision;
    op.precisionClass = precisionClass;

    op.synthesis = @synthesis;
    op.synthesisAdjoint = @synthesisAdjoint;
    op.forward = @forward;
    op.adjoint = @adjoint;
    op.forwardVector = @forwardVector;
    op.adjointVector = @adjointVector;
    op.normalVector = @normalVector;

    function X = synthesis(C)
        expectedSize = [nX nY nBasis];

        if ~isequal(size(C), expectedSize)
            error('buildGARSOSSTMOperators:SynthesisSize', ...
                'C must have size [%d %d %d].', ...
                nX, nY, nBasis);
        end

        if useSinglePrecision
            C = single(C);
        end

        X = sum(ST_maps .* ...
            reshape(C, [nX nY 1 nBasis]), 4);

        if useSinglePrecision && ~isa(X, 'single')
            X = single(X);
        end
    end

    function C = synthesisAdjoint(X)
        expectedSize = [nX nY nFrame];

        if ~isequal(size(X), expectedSize)
            error('buildGARSOSSTMOperators:SynthesisAdjointSize', ...
                'X must have size [%d %d %d].', ...
                nX, nY, nFrame);
        end

        if useSinglePrecision
            X = single(X);
        end

        C = sum(conjugateSTMaps .* ...
            reshape(X, [nX nY nFrame 1]), 3);

        if useSinglePrecision && ~isa(C, 'single')
            C = single(C);
        end
    end

    function y = forward(C)
        if useSinglePrecision
            C = single(C);
        end

        X = synthesis(C);
        y = radialForward(X);

        if useSinglePrecision && ~isa(y, 'single')
            y = single(y);
        end
    end

    function C = adjoint(y)
        if useSinglePrecision
            y = single(y);
        end

        X = radialAdjoint(y);
        C = synthesisAdjoint(X);

        if useSinglePrecision && ~isa(C, 'single')
            C = single(C);
        end
    end

    function yVector = forwardVector(cVector)
        if useSinglePrecision
            cVector = single(cVector);
        end

        C = reshape(cVector, [nX nY nBasis]);
        y = forward(C);
        yVector = y(:);
    end

    function cVector = adjointVector(yVector)
        if useSinglePrecision
            yVector = single(yVector);
        end

        y = reshape(yVector, ...
            [nReadout, nSpokes, nCoil, nFrame]);

        C = adjoint(y);
        cVector = C(:);
    end

    function outputVector = normalVector(inputVector)
        if useSinglePrecision
            inputVector = single(inputVector);
        end

        outputVector = adjointVector( ...
            forwardVector(inputVector));

        if useSinglePrecision && ~isa(outputVector, 'single')
            outputVector = single(outputVector);
        end
    end
end