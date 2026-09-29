function preconditioner = buildGARSOSSTMPreconditioner( ...
    stmOperator, varargin)
%buildGARSOSSTMPreconditioner Build an STM-domain diagonal preconditioner.
%
% Optional name-value arguments:
%   'Lambda'           Tikhonov parameter. Default: 0.
%   'NUFFTApproximation'
%                      'scalar' or 'image'. Default: 'scalar'.
%   'FloorFraction'    Minimum diagonal as a fraction of the maximum.
%                      Default: 1e-4.
%   'Normalize'        Normalize the NUFFT diagonal estimates. Default true.
%
% Returned fields:
%   diagonal           [Nx, Ny, nBasis] positive diagonal approximation
%   inverseDiagonal    reciprocal diagonal
%   apply              function handle implementing M \ r
%   nufftDiagonal      estimated per-frame NUFFT diagonal
%
% PCG usage:
%   Mfun = preconditioner.apply;
%   x = pcg(Afun, rhs, tol, maxIter, Mfun, [], x0);

    p = inputParser;
    p.FunctionName = mfilename;

    addParameter(p, 'Lambda', 0, ...
        @(x) isnumeric(x) && isscalar(x) && isreal(x) && x >= 0);

    addParameter(p, 'NUFFTApproximation', 'scalar', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));

    addParameter(p, 'FloorFraction', 1e-4, ...
        @(x) isnumeric(x) && isscalar(x) && ...
        isreal(x) && x > 0 && x < 1);

    addParameter(p, 'Normalize', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));

    parse(p, varargin{:});

    lambda = p.Results.Lambda;
    approximation = lower(string(p.Results.NUFFTApproximation));
    floorFraction = p.Results.FloorFraction;
    normalizeEstimate = logical(p.Results.Normalize);

    if ~ismember(approximation, ["scalar", "image"])
        error('buildGARSOSSTMPreconditioner:InvalidApproximation', ...
            'NUFFTApproximation must be ''scalar'' or ''image''.');
    end

    requiredFields = { ...
        'imageSize', 'nFrame', 'nBasis', 'ST_maps', ...
        'senseMaps', 'radialOperator'};

    if ~all(isfield(stmOperator, requiredFields))
        error('buildGARSOSSTMPreconditioner:InvalidOperator', ...
            'stmOperator is missing required fields.');
    end

    radialOperator = stmOperator.radialOperator;

    if ~isfield(radialOperator, 'plans')
        error('buildGARSOSSTMPreconditioner:MissingPlans', ...
            'The radial operator must expose its NUFFT plans.');
    end

    nx = stmOperator.imageSize(1);
    ny = stmOperator.imageSize(2);
    nFrame = stmOperator.nFrame;
    nBasis = stmOperator.nBasis;

    STMaps = stmOperator.ST_maps;
    senseMaps = stmOperator.senseMaps;

    if isfield(stmOperator, 'useSinglePrecision') && ...
            stmOperator.useSinglePrecision
        numericClass = 'single';
    else
        numericClass = class(STMaps);
    end

    coilPower = sum(abs(senseMaps).^2, 3);
    coilPower = cast(coilPower, numericClass);

    if approximation == "scalar"
        nufftDiagonal = zeros(1, 1, nFrame, numericClass);
    else
        nufftDiagonal = zeros(nx, ny, nFrame, numericClass);
    end

    % Estimate diag(F_t^H F_t) by applying the single-coil NUFFT normal
    % operator to an image of ones.
    onesImage = ones(nx, ny, numericClass);

    for t = 1:nFrame
        forwardValues = nufft(onesImage, radialOperator.plans(t).st);
        normalOnes = nufft_adj(forwardValues, ...
            radialOperator.plans(t).st);

        normalOnes = real(normalOnes);
        normalOnes = max(normalOnes, 0);
        normalOnes = cast(normalOnes, numericClass);

        if approximation == "scalar"
            positiveValues = normalOnes(normalOnes > 0);

            if isempty(positiveValues)
                frameScale = cast(1, numericClass);
            else
                frameScale = median(positiveValues(:));
            end

            nufftDiagonal(1, 1, t) = frameScale;
        else
            nufftDiagonal(:, :, t) = normalOnes;
        end
    end

    if normalizeEstimate
        positiveValues = nufftDiagonal(nufftDiagonal > 0);

        if ~isempty(positiveValues)
            scale = median(positiveValues(:));
            if scale > 0
                nufftDiagonal = nufftDiagonal ./ scale;
            end
        end
    end

    diagonal = zeros(nx, ny, nBasis, numericClass);

    for basisIndex = 1:nBasis
        basisPower = abs(STMaps(:, :, :, basisIndex)).^2;

        weightedBasisPower = basisPower .* nufftDiagonal;

        diagonal(:, :, basisIndex) = coilPower .* ...
            sum(weightedBasisPower, 3);
    end

    lambdaTyped = cast(lambda, numericClass);
    diagonal = real(diagonal) + lambdaTyped;

    maximumDiagonal = max(diagonal(:));

    if ~isfinite(maximumDiagonal) || maximumDiagonal <= 0
        error('buildGARSOSSTMPreconditioner:InvalidDiagonal', ...
            'The estimated preconditioner diagonal is not positive.');
    end

    floorValue = cast(floorFraction, numericClass) .* maximumDiagonal;
    diagonal = max(diagonal, floorValue);

    inverseDiagonal = 1 ./ diagonal;

    preconditioner = struct();
    preconditioner.diagonal = diagonal;
    preconditioner.inverseDiagonal = inverseDiagonal;
    preconditioner.nufftDiagonal = nufftDiagonal;
    preconditioner.coilPower = coilPower;
    preconditioner.lambda = lambda;
    preconditioner.approximation = char(approximation);
    preconditioner.floorFraction = floorFraction;
    preconditioner.apply = @applyPreconditioner;

    function outputVector = applyPreconditioner(inputVector)
        inputArray = reshape(inputVector, [nx, ny, nBasis]);
        outputArray = inverseDiagonal .* inputArray;
        outputVector = outputArray(:);
    end
end