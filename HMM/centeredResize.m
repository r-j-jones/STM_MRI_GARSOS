function newArray = centeredResize(oldArray,newArraySize)

if iscell(oldArray)
    
    newArray = cellfun(@(x)centeredResize(x,newArraySize), oldArray, 'UniformOutput', false);
    
else
    
    nDimensions = max(ndims(oldArray),numel(newArraySize));
    
    newArraySize = [newArraySize, ones(1, nDimensions - numel(newArraySize))];
    
    oldArraySize = [size(oldArray), ones(1, nDimensions - ndims(oldArray))];
    
    subArraySize = min(oldArraySize, newArraySize);
    
    oldArrayCenter = floor(oldArraySize/2) + 1;
    newArrayCenter = floor(newArraySize/2) + 1;
    subArrayCenter = floor(subArraySize/2) + 1;
    
    newArray = zeros(newArraySize,class(oldArray));
    
    oldStartIndices = 1 + (oldArrayCenter - 1) - (subArrayCenter - 1);
    oldEndIndices = oldStartIndices + (subArraySize - 1);
    
    newStartIndices = 1 + (newArrayCenter - 1) - (subArrayCenter - 1);
    newEndIndices = newStartIndices + (subArraySize - 1);
    
    oldIndices = arrayfun(@(x,y)x:y,oldStartIndices,oldEndIndices,'UniformOutput',false);
    newIndices = arrayfun(@(x,y)x:y,newStartIndices,newEndIndices,'UniformOutput',false);
    
    newArray(newIndices{:}) = oldArray(oldIndices{:});
    
end