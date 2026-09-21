function transformedSignal = centeredIFFT(untransformedSignal,nDimensions)

if iscell(untransformedSignal)
    
    transformedSignal = cellfun(@(x)centeredIFFT(x,nDimensions), untransformedSignal, 'UniformOutput', false);
    
else
    
    transformedSignal = zeros(size(untransformedSignal));
    
    nArrayDimensions = max(nDimensions, ndims(untransformedSignal));
    arraySize = cell(1, nArrayDimensions);
    [arraySize{:}] = size(untransformedSignal);
    arraySize = cell2mat(arraySize);
    nCoils = prod(arraySize)/prod(arraySize(1:nDimensions));
    arraySize = arraySize(1:nDimensions);
    
    for i = 1:nCoils
        indices = arrayfun(@(x,y)x:y, [ones(1,numel(arraySize)), i], [arraySize, i],'UniformOutput',false);
        transformedSignal(indices{:}) = fftshift(ifftn(ifftshift(untransformedSignal(indices{:}))));
    end
    
end



