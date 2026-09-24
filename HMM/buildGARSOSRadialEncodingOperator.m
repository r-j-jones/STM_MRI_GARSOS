function op = buildGARSOSRadialEncodingOperator(kxy, imageSize, senseMaps, varargin)
%buildGARSOSRADIALENCODINGOPERATOR Build matched multicoil radial operators.
%
% Forward:
%   F(X)(:, :, c, t) = NUFFT_t(S_c .* X(:, :, t))
% Adjoint:
%   F^H(y)(:, :, t) = sum_c conj(S_c) .* NUFFT_t^H(y(:, :, c, t))
%
% No density compensation, empirical scaling, coil normalization, or
% automatic coil combination is applied.

    validateattributes(imageSize, {'numeric'}, {'vector', 'numel', 2, 'integer', 'positive'});
    validateattributes(kxy, {'numeric'}, {'3d', 'nonempty'});
    validateattributes(senseMaps, {'numeric'}, {'3d', 'nonempty'});
    if size(senseMaps, 1) ~= imageSize(1) || ...
            size(senseMaps, 2) ~= imageSize(2)
        error('buildGARSOSRadialEncodingOperator:SenseMapSize', ...
            'senseMaps must have spatial size imageSize.');
    end

    nCoil = size(senseMaps, 3);
    nFrame = size(kxy, 3);
    plans = buildGARSOSRadialNUFFTPlans(kxy, imageSize, varargin{:});
    nReadout = size(kxy, 1);
    nSpokes = size(kxy, 2);

    op = struct();
    op.imageSize = imageSize(:).';
    op.nReadout = nReadout;
    op.nSpokes = nSpokes;
    op.nCoil = nCoil;
    op.nFrame = nFrame;
    op.senseMaps = senseMaps;
    op.plans = plans;
    op.forward = @forward;
    op.adjoint = @adjoint;
    op.forwardVector = @(x) forward(reshape(x, [imageSize, nFrame]));
    op.adjointVector = @(y) reshape(adjoint(reshape(y, [nReadout, nSpokes, nCoil, nFrame])), [], 1);

    function y = forward(x)
        validateattributes(x, {'numeric'}, {'size', [imageSize, nFrame]});
        y = zeros(nReadout, nSpokes, nCoil, nFrame, 'like', x);
        for t = 1:nFrame
            for c = 1:nCoil
                imageCoil = x(:, :, t) .* senseMaps(:, :, c);
                y(:, :, c, t) = reshape( ...
                    nufft(imageCoil, plans(t).st), [nReadout, nSpokes]);
            end
        end
    end

    function x = adjoint(y)
        validateattributes(y, {'numeric'}, ...
            {'size', [nReadout, nSpokes, nCoil, nFrame]});
        x = zeros(imageSize(1), imageSize(2), nFrame, 'like', y);
        for t = 1:nFrame
            imageAdjoint = zeros(imageSize, 'like', y);
            for c = 1:nCoil
                coilAdjoint = nufft_adj( ...
                    reshape(y(:, :, c, t), [], 1), plans(t).st);
                imageAdjoint = imageAdjoint + ...
                    conj(senseMaps(:, :, c)) .* coilAdjoint;
            end
            x(:, :, t) = imageAdjoint;
        end
    end
end
