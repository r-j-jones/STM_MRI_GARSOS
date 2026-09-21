function kernelSamples = generateKaiserBesselKernel(imageWidth, gridWidth, kernelWidth, kernelSamplingDensity)

assert(mod(kernelWidth*kernelSamplingDensity/2, 1) == 0, 'Kernel must be sampled with an 2*integer number of samples');

nSamples = kernelWidth*kernelSamplingDensity/2 + 1;

kernelRadii = linspace(0,kernelWidth/2/gridWidth,nSamples);

alpha = gridWidth/imageWidth;

beta = pi*sqrt((kernelWidth/alpha)^2*(alpha - 0.5)^2 - 0.8);

kernelSamples = (gridWidth/kernelWidth)*besseli(0, beta*sqrt(1 - (2*gridWidth*kernelRadii/kernelWidth).^2));

kernelSamples = kernelSamples/(trapz([fliplr(kernelSamples(2:end)) kernelSamples])/kernelSamplingDensity);