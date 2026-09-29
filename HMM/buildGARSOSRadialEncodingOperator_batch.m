function op = buildGARSOSRadialEncodingOperator( ...
    kxy, imageSize, senseMaps, varargin)
%buildGARSOSRADIALENCODINGOPERATOR Build matched multicoil radial operators.
%
% Forward:
%   F(X)(:,:,c,t) = NUFFT_t(S_c .* X(:,:,t))
%
% Adjoint:
%   F^H(Y)(:,:,t) =
%       sum_c conj(S_c) .* NUFFT_t^H(Y(:,:,c,t))
%
% Coils are processed as a batch for each temporal frame.

    validateattributes(imageSize, {'numeric'}, ...
        {'vector', 'numel', 2, 'integer', 'positive'});

    validateattributes(kxy, {'numeric'}, {'nonempty'});
    validateattributes(senseMaps, {'numeric'}, {'nonempty'});

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

    nReadout = size(kxy, 1);
    nSpokes = size(kxy, 2);
    nFrame = size(kxy, 3);
    nCoil = size(senseMaps, 3);
    nSamplesPerFrame = nReadout * nSpokes;

    plans = buildGARSOSRadialNUFFTPlans( ...
        kxy, imageSize, varargin{:});

    % Avoid recomputing the conjugate during every adjoint operation.
    conjugateSenseMaps = conj(senseMaps);

    op = struct();
    op.imageSize = imageSize;
    op.nReadout = nReadout;
    op.nSpokes = nSpokes;
    op.nCoil = nCoil;
    op.nFrame = nFrame;
    op.senseMaps = senseMaps;
    op.plans = plans;

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

        y = zeros( ...
            nReadout, nSpokes, nCoil, nFrame, 'like', x);

        for t = 1:nFrame
            % Implicitly expand the frame over all coil maps.
            % Size: [Nx, Ny, nCoil]
            coilImages = x(:, :, t) .* senseMaps;

            % Batched NUFFT.
            % Input:  [Nx, Ny, nCoil]
            % Output: [nReadout*nSpokes, nCoil]
            frameSamples = nufft(coilImages, plans(t).st);

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

        x = zeros(imageSize(1), imageSize(2), nFrame, 'like', y);

        for t = 1:nFrame
            % One column for each coil.
            % Size: [nReadout*nSpokes, nCoil]
            frameSamples = reshape( ...
                y(:, :, :, t), nSamplesPerFrame, nCoil);

            % Batched adjoint NUFFT.
            % Output: [Nx, Ny, nCoil]
            coilAdjoint = nufft_adj( ...
                frameSamples, plans(t).st);

            % Sensitivity-weighted coil combination.
            x(:, :, t) = sum( ...
                conjugateSenseMaps .* coilAdjoint, 3);
        end
    end
end