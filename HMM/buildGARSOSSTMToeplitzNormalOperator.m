function op = buildGARSOSSTMToeplitzNormalOperator( ...
    stmOperator, radialToeplitzNormal)
%buildGARSOSSTMToeplitzNormalOperator Compose STM with radial Toeplitz Gram.
%
% Applies:
%   T_STM^H * A_radial^H A_radial * T_STM

    imageSize = stmOperator.imageSize;
    nBasis = stmOperator.nBasis;

    if ~isequal(radialToeplitzNormal.imageSize, imageSize) || ...
            radialToeplitzNormal.nFrame ~= stmOperator.nFrame
        error('buildGARSOSSTMToeplitzNormalOperator:DimensionMismatch', ...
            'STM and radial Toeplitz dimensions disagree.');
    end

    synthesis = stmOperator.synthesis;
    synthesisAdjoint = stmOperator.synthesisAdjoint;
    radialNormal = radialToeplitzNormal.apply;

    op = struct();
    op.imageSize = imageSize;
    op.nBasis = nBasis;
    op.apply = @apply;
    op.applyVector = @applyVector;

    function output = apply(input)
        dynamicImages = synthesis(input);
        dynamicNormal = radialNormal(dynamicImages);
        output = synthesisAdjoint(dynamicNormal);
    end

    function outputVector = applyVector(inputVector)
        input = reshape(inputVector, [imageSize, nBasis]);
        output = apply(input);
        outputVector = output(:);
    end
end