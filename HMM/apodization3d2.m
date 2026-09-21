function [gridApodizationFunction, imageApodizationFunction] = apodization3d2(gridSize,imageSize,kernelSize)

kernelSamplingDensity = 200;

nDimensions = numel(gridSize);

gridApodizationFunction = cell(nDimensions,1);
imageApodizationFunction = cell(nDimensions,1);

for iDimension = 1:nDimensions
    
    gridWidth = gridSize(iDimension);
    imageWidth = imageSize(iDimension);
    kernelWidth = kernelSize(iDimension);    
    
    kernelSamples = generateKaiserBesselKernel(imageWidth, gridWidth, kernelWidth, kernelSamplingDensity);
    
    kernelSamples = centeredResize([fliplr(kernelSamples(2:end)) kernelSamples],[1 kernelSamplingDensity*gridWidth]);
    
    kernelSamples = fftshift(fftn(ifftshift(kernelSamples)));
    
    cx = centeredResize(kernelSamples, [1 gridWidth]);
    hx = sinc2((-floor(gridWidth/2):floor((gridWidth - 1)/2))/(kernelSamplingDensity*gridWidth));
    
    cx = cx.*hx;
    
    ccx = shiftAndResize(cx, [], [1 imageWidth], 'mm', 'mm');
    
    gridApodizationFunction{iDimension} = real(shiftdim(cx(:)/kernelSamplingDensity, 1 - iDimension));
    imageApodizationFunction{iDimension} = real(shiftdim(ccx(:)/kernelSamplingDensity, 1 - iDimension));
    
end

% cx = apodization1D(kernelSamplingDensity,gridWidth,imageWidth,kernelWidth);

%Y = reshape(reshape(cx'*cx,[imageWidth*imageWidth 1])*cx,[imageWidth imageWidth imageWidth])/(kernelSamplingDensity^3);