function Z2 = decorrelateCoils2(Z, noiseWhiteningTransform)

%Z2 = tprod(Z, [1 2 3 -5], noiseWhiteningTransform, [-5 4]);

% pg = Progressor('Decorrelating coils.');
% 
% Z2 = zeros([size(Z, 1), size(Z, 2), size(Z, 3), size(noiseWhiteningTransform, 2)], 'like', Z);
% 
% nElements = size(Z, 4);
% 
% assert(nElements == size(noiseWhiteningTransform, 1));
% 
% for iElement = 1:nElements
% 
%     Z2 = Z2 + bsxfun(@times, Z(:, :, :, iElement), permute(noiseWhiteningTransform(iElement, :), [3 4 1 2]));
% 
%     pg.setProgress(iElement/nElements);
%     
% end

arraySize = cell(1,4);

[arraySize{:}] = size(Z);

arraySize = cell2mat(arraySize);

assert(arraySize(4) == size(noiseWhiteningTransform,1));

Z2 = reshape(reshape(Z, [prod(arraySize(1:3)), arraySize(4)])*noiseWhiteningTransform, [arraySize(1:3), size(noiseWhiteningTransform,2)]);


