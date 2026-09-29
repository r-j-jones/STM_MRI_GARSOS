function op = buildGARSOSRadialToeplitzNormalOperator( ...
    radialOperator, toeplitzGram)
%buildGARSOSRadialToeplitzNormalOperator Build multicoil Toeplitz normal.
%
% Applies, for each frame:
%   sum_c conj(S_c) .* T_t(S_c .* x_t)
%
% where T_t approximates or represents F_t^H F_t.

    imageSize = radialOperator.imageSize;
    nFrame = radialOperator.nFrame;
    nCoil = radialOperator.nCoil;
    senseMaps = radialOperator.senseMaps;
    conjugateSenseMaps = conj(senseMaps);

    if toeplitzGram.nFrame ~= nFrame || ...
            ~isequal(toeplitzGram.imageSize, imageSize)
        error('buildGARSOSRadialToeplitzNormalOperator:DimensionMismatch', ...
            'Toeplitz Gram and radial operator dimensions disagree.');
    end

    op = struct();
    op.imageSize = imageSize;
    op.nFrame = nFrame;
    op.nCoil = nCoil;
    op.apply = @apply;
    op.applyFrame = @applyFrame;

    function output = apply(input)
        if ~isequal(size(input), [imageSize, nFrame])
            error('buildGARSOSRadialToeplitzNormalOperator:InputSize', ...
                'Input must have size [%d %d %d].', ...
                imageSize(1), imageSize(2), nFrame);
        end

        output = zeros(size(input), 'like', input);

        for t = 1:nFrame
            output(:, :, t) = applyFrame(input(:, :, t), t);
        end
    end

    function outputFrame = applyFrame(inputFrame, t)
        % Nx x Ny x nCoil
        coilImages = inputFrame .* senseMaps;

        % Apply the frame Gram to all coils as a batch.
        coilNormal = toeplitzGram.applyFrame(coilImages, t);

        % Some IRT matrix classes may return double.
        if isa(inputFrame, 'single') && ~isa(coilNormal, 'single')
            coilNormal = single(coilNormal);
        end

        outputFrame = sum( ...
            conjugateSenseMaps .* coilNormal, 3);
    end
end