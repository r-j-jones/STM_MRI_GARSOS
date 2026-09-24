classdef GriddingReconstructur
    
    properties (Access  = public)
        
        apodizationFunction = [];
        II = [];
        DK = [];
        nSamples = [];
        imageSize = [];
        gridSize = [];
        kernelSize = [];
        nDimensions = [];
        
    end
    
    methods
        
        function obj = GriddingReconstructur(gridSize, imageSize, kernelSize, sampleCoordinates, coilCombinationProfiles)
            
            nDimensions = numel(gridSize);
            
            nSamples = size(sampleCoordinates, 1);
            
            kernelSamplingDensity = 200;
            
            kernelSamples = arrayfun(@(imageWidth, gridWidth, kernelWidth)generateKaiserBesselKernel(imageWidth, gridWidth, kernelWidth, kernelSamplingDensity), imageSize, gridSize, kernelSize, 'UniformOutput', false);
            
            [~, apodizationFunctions] = apodization3d2(gridSize, imageSize, kernelSize);
            
            apodizationFunction = coilCombinationProfiles*prod(gridSize);
            
            for iDimension = 1:nDimensions
                
                apodizationFunction = apodizationFunction./apodizationFunctions{iDimension};
                
            end
            
            apodizationFunction = reshape(apodizationFunction, prod(imageSize), []);
            
            kernelSamples2 = cellfun(@(x)[fliplr(x(1, 2:end)), x], kernelSamples, 'UniformOutput', false);
            
            kKernel = cellfun(@(x)fftlattice(1/kernelSamplingDensity, 0, numel(x)), kernelSamples2, 'UniformOutput', false);
            
            p = arrayfun(@(x)-floor(x/2):floor(x/2), kernelSize, 'UniformOutput', false);
            
            pp = fftlattice(eye(nDimensions), zeros(nDimensions, 1), kernelSize);
            
            kGrid = sampleCoordinates.*gridSize;
            
            ii = round(kGrid);
            
            dk = ii - kGrid;
            
            gridStride = [1 cumprod(gridSize(1, 1:end-1))];
            
            II = mod(reshape(permute(ii, [1 3 2]) + permute(pp, [3 2 1]), [nSamples*prod(kernelSize), nDimensions]), gridSize)*gridStride.' + 1;
            
            for iDimension = 1:nDimensions
                
                dimOrder = [1 circshift([2 3 4], iDimension - 1)];
                
                if iDimension == 1
                    
                    DK = permute(interp1(kKernel{iDimension}, kernelSamples2{iDimension}, dk(:, iDimension) + p{iDimension}), dimOrder);
                    
                else
                    
                    DK = DK.*permute(interp1(kKernel{iDimension}, kernelSamples2{iDimension}, dk(:, iDimension) + p{iDimension}), dimOrder);
                    
                end
                
            end
            
            DK = reshape(DK, nSamples, prod(kernelSize));
            
            obj.II = II;
            obj.DK = DK;
            obj.nSamples = nSamples;
            obj.gridSize = gridSize;
            obj.imageSize = imageSize;
            obj.kernelSize = kernelSize;
            obj.nDimensions = nDimensions;
            obj.apodizationFunction = apodizationFunction;
            
        end
        
        function imageData = grid(obj, sampledData, outputIsGpuArray)
            nCurrentSamples = size(sampledData, 1);
            assert(obj.nSamples == nCurrentSamples);
            nChannels = size(sampledData, 2);
            for iChannel = 1:1:nChannels
                DKi = reshape(obj.DK.*sampledData(:, iChannel), [obj.nSamples*prod(obj.kernelSize), 1]);
                A = accumarray(obj.II, DKi, [prod(obj.gridSize), 1]);
                A = reshape(A, obj.gridSize);
                A = ifftn(A);
                A = shiftAndResize(A, [], obj.imageSize, repmat('e', [1, obj.nDimensions]), repmat('m', [1, obj.nDimensions]));
                A = A(:);
                A = A.*obj.apodizationFunction(:,iChannel);
                if iChannel == 1
                    imageData = A;
                else
                   imageData = imageData + A;
                end                
            end
            imageData = reshape(imageData, obj.imageSize);
            if ~outputIsGpuArray
                imageData = gather(imageData);
            end       
        end


        % Implement method to return gridded coil-resolved k-space
        function kData = gridNoCombine(obj, sampledData, outputIsGpuArray)
            nCurrentSamples = size(sampledData, 1);
            assert(obj.nSamples == nCurrentSamples);
            nChannels = size(sampledData, 2);
            for iChannel = 1:1:nChannels
                DKi = reshape(obj.DK.*sampledData(:, iChannel), [obj.nSamples*prod(obj.kernelSize), 1]);
                A = accumarray(obj.II, DKi, [prod(obj.gridSize), 1]);
                A = reshape(A, obj.gridSize);
                A = ifftn(A);
                A = shiftAndResize(A, [], obj.imageSize, repmat('e', [1, obj.nDimensions]), repmat('m', [1, obj.nDimensions]));
                A = A(:);
                A = A.*obj.apodizationFunction(:,iChannel);
                if iChannel == 1
                    kData = A;
                else
                   kData = kData + A;
                end                
            end
            kData = reshape(kData, obj.imageSize);
            if ~outputIsGpuArray
                kData = gather(kData);
            end       
        end


        
    end
    
end