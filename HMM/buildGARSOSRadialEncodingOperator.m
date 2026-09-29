function op = buildGARSOSRadialEncodingOperator( ...
    kxy, imageSize, senseMaps, varargin)
%buildGARSOSRadialEncodingOperator Build matched multicoil radial operators.
%
% Forward:
%   F(X)(:,:,c,t) = NUFFT_t(S_c .* X(:,:,t))
%
% Adjoint:
%   F^H(Y)(:,:,t) =
%       sum_c conj(S_c) .* NUFFT_t^H(Y(:,:,c,t))
%
% Coils are processed as a batch for each temporal frame.
%
% Optional name-value argument:
%   'UseSinglePrecision'  Convert trajectory, sensitivity maps, operator
%                         inputs, and operator outputs to single precision.
%                         Default: false.
%
% All remaining name-value arguments are forwarded to
% buildGARSOSRadialNUFFTPlans.
%
% Example:
%   op = buildGARSOSRadialEncodingOperator( ...
%       kxy, imageSize, senseMaps, ...
%       'UseSinglePrecision', true);

    validateattributes(imageSize, {'numeric'}, ...
        {'vector', 'numel', 2, 'integer', 'positive'});

    validateattributes(kxy, {'numeric'}, {'nonempty'});
    validateattributes(senseMaps, {'numeric'}, {'nonempty'});

    % Parse this function's option while preserving all other options for
    % buildGARSOSRadialNUFFTPlans.
    [useSinglePrecision, planArguments] = ...
        extractUseSinglePrecision(varargin);

    imageSize = imageSize(:).';

    if ndims(kxy) ~= 3
        error('buildGARSOSRadialEncodingOperator:KxyDimensions', ...
            'kxy must have size [nReadout, nSpokes, nFrame].');
    end

    if ndims(senseMaps) ~= 3
        error('buildGARSOSRadialEncodingOperator:SenseDimensions', ...
            'senseMaps must have size [Nx, Ny, nCoil].');
    end

    if size(senseMaps, 1) ~= imageSize(1) || ...
            size(senseMaps, 2) ~= imageSize(2)
        error('buildGARSOSRadialEncodingOperator:SenseMapSize', ...
            'senseMaps must have spatial size imageSize.');
    end

    if useSinglePrecision
        kxy = single(kxy);
        senseMaps = single(senseMaps);
        precisionClass = 'single';
    else
        precisionClass = class(senseMaps);
    end

    nReadout = size(kxy, 1);
    nSpokes = size(kxy, 2);
    nFrame = size(kxy, 3);
    nCoil = size(senseMaps, 3);
    nSamplesPerFrame = nReadout * nSpokes;

    plans = buildGARSOSRadialNUFFTPlans( ...
        kxy, imageSize, planArguments{:});

    % Convert public plan arrays where possible. The internal MIRT st
    % structure may still contain doubles depending on the implementation.
    if useSinglePrecision
        plans = convertPlanNumericFieldsToSingle(plans);
    end

    conjugateSenseMaps = conj(senseMaps);

    op = struct();
    op.imageSize = imageSize;
    op.nReadout = nReadout;
    op.nSpokes = nSpokes;
    op.nCoil = nCoil;
    op.nFrame = nFrame;
    op.senseMaps = senseMaps;
    op.plans = plans;
    op.useSinglePrecision = useSinglePrecision;
    op.precisionClass = precisionClass;

    op.forward = @forward;
    op.adjoint = @adjoint;

    op.forwardVector = @(x) ...
        forward(reshape(x, [imageSize, nFrame]));

    op.adjointVector = @(y) reshape( ...
        adjoint(reshape(y, ...
        [nReadout, nSpokes, nCoil, nFrame])), [], 1);

    function y = forward(x)
        expectedSize = [imageSize, nFrame];

        if ~isequal(size(x), expectedSize)
            error('buildGARSOSRadialEncodingOperator:ForwardSize', ...
                'x must have size [%d %d %d].', ...
                imageSize(1), imageSize(2), nFrame);
        end

        if useSinglePrecision
            x = single(x);
        end

        y = zeros( ...
            nReadout, nSpokes, nCoil, nFrame, ...
            'like', x);

        for t = 1:nFrame
            % Size: [Nx, Ny, nCoil]
            coilImages = x(:, :, t) .* senseMaps;

            % Output: [nReadout*nSpokes, nCoil]
            frameSamples = nufft(coilImages, plans(t).st);

            % Some NUFFT implementations promote single inputs to double.
            if useSinglePrecision && ~isa(frameSamples, 'single')
                frameSamples = single(frameSamples);
            end

            y(:, :, :, t) = reshape( ...
                frameSamples, nReadout, nSpokes, nCoil);
        end
    end

    function x = adjoint(y)
        expectedSize = [nReadout, nSpokes, nCoil, nFrame];

        if ~isequal(size(y), expectedSize)
            error('buildGARSOSRadialEncodingOperator:AdjointSize', ...
                'y must have size [%d %d %d %d].', ...
                nReadout, nSpokes, nCoil, nFrame);
        end

        if useSinglePrecision
            y = single(y);
        end

        x = zeros(imageSize(1), imageSize(2), nFrame, 'like', y);

        for t = 1:nFrame
            frameSamples = reshape( ...
                y(:, :, :, t), nSamplesPerFrame, nCoil);

            coilAdjoint = nufft_adj( ...
                frameSamples, plans(t).st);

            if useSinglePrecision && ~isa(coilAdjoint, 'single')
                coilAdjoint = single(coilAdjoint);
            end

            x(:, :, t) = sum( ...
                conjugateSenseMaps .* coilAdjoint, 3);
        end
    end
end


function [useSinglePrecision, remainingArguments] = ...
    extractUseSinglePrecision(arguments)
% Extract UseSinglePrecision without consuming NUFFT-plan options.

    useSinglePrecision = false;
    remainingArguments = {};

    if mod(numel(arguments), 2) ~= 0
        error('buildGARSOSRadialEncodingOperator:InvalidOptions', ...
            'Optional inputs must be supplied as name-value pairs.');
    end

    optionFound = false;

    for index = 1:2:numel(arguments)
        optionName = arguments{index};
        optionValue = arguments{index + 1};

        if ischar(optionName) || ...
                (isstring(optionName) && isscalar(optionName))

            if strcmpi(string(optionName), "UseSinglePrecision")
                if optionFound
                    error(['buildGARSOSRadialEncodingOperator:' ...
                        'DuplicatePrecisionOption'], ...
                        'UseSinglePrecision was specified more than once.');
                end

                validateattributes(optionValue, {'logical', 'numeric'}, ...
                    {'scalar'}, mfilename, 'UseSinglePrecision');

                if isnumeric(optionValue) && ...
                        ~(optionValue == 0 || optionValue == 1)
                    error(['buildGARSOSRadialEncodingOperator:' ...
                        'InvalidPrecisionOption'], ...
                        'UseSinglePrecision must be true or false.');
                end

                useSinglePrecision = logical(optionValue);
                optionFound = true;
                continue;
            end
        end

        remainingArguments(end+1:end+2) = ... %#ok<AGROW>
            {optionName, optionValue};
    end
end


function plans = convertPlanNumericFieldsToSingle(plans)
% Convert ordinary numeric plan fields to single recursively.
%
% MIRT's internal st structure can rely on specific field types. To avoid
% altering lookup-table indices or sparse matrices incorrectly, only known
% floating-point trajectory metadata fields are converted here. The NUFFT
% builder should ideally create an st structure suitable for its own
% implementation.

    fieldsToConvert = {'normalizedTrajectory', 'omega'};

    for planIndex = 1:numel(plans)
        for fieldIndex = 1:numel(fieldsToConvert)
            fieldName = fieldsToConvert{fieldIndex};

            if isfield(plans(planIndex), fieldName) && ...
                    isnumeric(plans(planIndex).(fieldName))
                plans(planIndex).(fieldName) = ...
                    single(plans(planIndex).(fieldName));
            end
        end
    end
end