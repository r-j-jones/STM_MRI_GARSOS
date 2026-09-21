function Y = shiftAndResize(X, copySize, outputSize, inputZeroLocations, outputZeroLocations)

nDimensions = max([ndims(X), numel(copySize), numel(outputSize)]);

inputSize = getArraySize(X, nDimensions);

if isempty(outputSize)
    
    outputSize = inputSize;
    
else
    
    outputSize = [reshape(outputSize,1,[]), ones(1, nDimensions - numel(outputSize))];
    
end

if isempty(copySize)    
    
    copySize = inf(1, nDimensions);    
    
else
    
    copySize = [copySize, ones(1, nDimensions - numel(copySize))];    
    
end

if numel(inputZeroLocations) < nDimensions
    
    switch class(inputZeroLocations)
        
        case 'double'
            
            inputZeroLocations = [inputZeroLocations, ones(1, nDimensions - numel(inputZeroLocations))];
            
        case 'char'
            
            inputZeroLocations = [inputZeroLocations, repmat('e', [1, nDimensions - numel(inputZeroLocations)])];
            
        otherwise
            
            error('Input zeros location symbol must be double or char.');
            
    end
    
end

if numel(outputZeroLocations) < nDimensions
    
    switch class(outputZeroLocations)
        
        case 'double'
            
            outputZeroLocations = [outputZeroLocations, ones(1, nDimensions - numel(outputZeroLocations))];
            
        case 'char'
            
            outputZeroLocations = [outputZeroLocations, repmat('e', [1, nDimensions - numel(outputZeroLocations)])];
            
        otherwise
            
            error('Output zeros location symbol must be double or char.');
            
    end
    
end

copySize = min([copySize; inputSize; outputSize], [], 1);

iiInputs = cell(1, nDimensions);
iiOutputs = cell(1, nDimensions);

for iDimension = 1:nDimensions
    
    inputWidth = inputSize(iDimension);
    cropWidth = copySize(iDimension);
    outputWidth = outputSize(iDimension);
    inputZeroLocation = inputZeroLocations(iDimension);
    outputZeroLocation = outputZeroLocations(iDimension);
    
    if ischar(inputZeroLocation)
        
        switch inputZeroLocation
            
            case 'e'
                
                inputZeroIndex = 1;
                
            case 'm'
                
                inputZeroIndex = floor(inputWidth/2) + 1;
                
            otherwise
                
                error('Unknown zero location symbol "%s". Use "m" for the center (middle) of the array and "e" for the edge.', inputZeroLocation);
                
        end
        
    else
        
        inputZeroIndex = inputZeroLocation;
        
    end
    
    if ischar(outputZeroLocation)
        
        switch outputZeroLocation
            
            case 'e'
                
                outputZeroIndex = 1;
                
            case 'm'
                
                outputZeroIndex = floor(outputWidth/2) + 1;
                
            otherwise
                
                error('Unknown zero location symbol "%s". Use "m" for the center (middle) of the array and "e" for the edge.', outputZeroLocation);
                
        end
        
    else
        
        outputZeroIndex = outputZeroLocation;
        
    end
    
    iiCrop = (1:cropWidth) - (floor(cropWidth/2) + 1);
    
    iiInput = mod((iiCrop - 0) + (inputZeroIndex - 1), inputWidth) + 1;
    
    iiOutput = mod((iiCrop - 0) + (outputZeroIndex - 1), outputWidth) + 1;
    
    [iiOutput, iiSorted] = sort(iiOutput, 'ascend');
    
    iiInput = iiInput(iiSorted);
    
    iiInputs{1, iDimension} = iiInput;
    iiOutputs{1, iDimension} = iiOutput;
    
end

if all(outputSize <= inputSize) && all(outputSize == copySize)
    
    Y = X(iiInputs{:});
    
else
    
    Y = zeros(outputSize,'like',X);
    
    Y(iiOutputs{:}) = X(iiInputs{:});
    
end











