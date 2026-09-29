function toeplitz = buildGARSOSToeplitzGramOperators( ...
    radialOperator, varargin)
%buildGARSOSToeplitzGramOperators Build framewise IRT Toeplitz Gram operators.
%
% Requires IRT functions:
%   Gnufft
%   build_gram
%
% Optional parameters:
%   'Jd'       NUFFT interpolation neighborhood. Default: [6 6]
%   'Kd'       Oversampled FFT size. Default: 2 * imageSize
%   'NShift'   Image-domain shift. Default: imageSize / 2
%   'Kernel'   NUFFT kernel option. Default: 'minmax:kb'
%   'Table'    Lookup table size. Default: 2^10
%   'Verbose'  Print progress. Default: true

    p = inputParser;
    p.FunctionName = mfilename;

    imageSize = radialOperator.imageSize;

    addParameter(p, 'Jd', [6 6], ...
        @(x) isnumeric(x) && numel(x) == 2);

    addParameter(p, 'Kd', 2 * imageSize, ...
        @(x) isnumeric(x) && numel(x) == 2);

    addParameter(p, 'NShift', imageSize / 2, ...
        @(x) isnumeric(x) && numel(x) == 2);

    addParameter(p, 'Kernel', 'minmax:kb', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));

    addParameter(p, 'Table', 2^10, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);

    addParameter(p, 'Verbose', true, ...
        @(x) (islogical(x) || isnumeric(x)) && isscalar(x));

    parse(p, varargin{:});

    Jd = double(p.Results.Jd(:).');
    Kd = double(p.Results.Kd(:).');
    nShift = double(p.Results.NShift(:).');
    kernel = char(p.Results.Kernel);
    tableSize = p.Results.Table;
    verbose = logical(p.Results.Verbose);

    requiredFields = {'imageSize', 'nFrame', 'plans'};

    if ~all(isfield(radialOperator, requiredFields))
        error('buildGARSOSToeplitzGramOperators:InvalidOperator', ...
            'radialOperator is missing required fields.');
    end

    if exist('Gnufft', 'file') == 0
        error('buildGARSOSToeplitzGramOperators:MissingGnufft', ...
            'Gnufft is not available on the MATLAB path.');
    end

    if exist('build_gram', 'file') == 0
        error('buildGARSOSToeplitzGramOperators:MissingBuildGram', ...
            'build_gram is not available on the MATLAB path.');
    end

    nFrame = radialOperator.nFrame;
    fullMask = true(imageSize);

    G = cell(1, nFrame);
    gram = cell(1, nFrame);

    nufftArguments = { ...
        imageSize, Jd, Kd, nShift, ...
        'table', tableSize, kernel};

    for t = 1:nFrame
        if verbose && (t == 1 || mod(t, 10) == 0 || t == nFrame)
            fprintf('Building Toeplitz Gram operator %d/%d...\n', ...
                t, nFrame);
        end

        omega = double(radialOperator.plans(t).omega);

        G{t} = Gnufft(fullMask, ...
            [{omega}, nufftArguments(:)']);

        gram{t} = build_gram(G{t}, 1);
    end

    toeplitz = struct();
    toeplitz.G = G;
    toeplitz.gram = gram;
    toeplitz.imageSize = imageSize;
    toeplitz.nFrame = nFrame;
    toeplitz.Jd = Jd;
    toeplitz.Kd = Kd;
    toeplitz.nShift = nShift;
    toeplitz.kernel = kernel;
    toeplitz.tableSize = tableSize;
    toeplitz.applyFrame = @applyFrame;

    function output = applyFrame(input, t)
        if ~isequal(size(input, 1), imageSize(1)) || ...
                ~isequal(size(input, 2), imageSize(2))
            error('buildGARSOSToeplitzGramOperators:InputSize', ...
                'Input spatial dimensions must equal imageSize.');
        end

        trailingSize = size(input);
        trailingSize = trailingSize(3:end);
        nBatch = prod(trailingSize);

        inputMatrix = reshape(input, prod(imageSize), nBatch);

        % A full true mask means the Gnufft domain has prod(imageSize)
        % elements. Tm can generally process multiple columns.
        outputMatrix = gram{t} * inputMatrix;

        output = reshape(outputMatrix, ...
            [imageSize, trailingSize]);
    end
end