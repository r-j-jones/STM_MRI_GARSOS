function [sliceData, diagnostics] = reconstructGARSOSRadialZSlice( ...
    sampledData, sampleCoordinates, nReadout, nPartition, nSpokes, ...
    coordinateTransform, imageSizeZ, gridSizeZ, kernelSizeZ, zIndex, varargin)
%reconstructGARSOSRADIALZSLICE Reconstruct radial data along the z axis.
%
% The input samples are restored as [readout, partition, spoke, coil].
% Only the transformed kz coordinate is gridded. Readout, spoke, coil,
% and in-plane trajectory ordering are preserved.

    p = inputParser;
    addParameter(p, 'TranslationVector', zeros(3, 1), ...
        @(x) isnumeric(x) && numel(x) == 3 && all(isfinite(x(:))));
    addParameter(p, 'KernelSamplingDensity', 200, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'KzTolerance', 1e-10, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(p, 'ReturnAllZ', false, ...
        @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    validateattributes(sampledData, {'numeric'}, {'2d', 'nonempty'});
    validateattributes(sampleCoordinates, {'numeric'}, {'2d', 'ncols', 3});
    validateattributes(nReadout, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(nPartition, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(nSpokes, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(coordinateTransform, {'numeric'}, ...
        {'size', [3 3], 'finite'});
    validateattributes(imageSizeZ, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(gridSizeZ, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(kernelSizeZ, {'numeric'}, ...
        {'scalar', 'integer', 'positive', 'odd'});
    validateattributes(zIndex, {'numeric'}, ...
        {'scalar', 'integer', 'positive', '<=', imageSizeZ});

    if size(sampledData, 1) ~= size(sampleCoordinates, 1)
        error('reconstructGARSOSRadialZSlice:SizeMismatch', ...
            'The sample and trajectory row counts must match.');
    end
    if any(~isfinite(sampledData(:))) || any(~isfinite(sampleCoordinates(:)))
        error('reconstructGARSOSRadialZSlice:NonFiniteInput', ...
            'Samples and coordinates must be finite.');
    end

    nSamples = size(sampledData, 1);
    nCoil = size(sampledData, 2);
    if nSamples ~= nReadout * nPartition * nSpokes
        error('reconstructGARSOSRadialZSlice:InvalidDimensions', ...
            'The flattened sample count is not compatible with the inferred dimensions.');
    end

    coordinates4 = reshape(sampleCoordinates, ...
        [nReadout, nPartition, nSpokes, 3]);
    data4 = reshape(sampledData, ...
        [nReadout, nPartition, nSpokes, nCoil]);

    translationVector = p.Results.TranslationVector(:);
    phase = exp(1i * 2 * pi * reshape(sum( ...
        coordinates4 .* reshape(translationVector, [1 1 1 3]), 4), ...
        [nReadout, nPartition, nSpokes]));
    data4 = data4 .* reshape(phase, [nReadout, nPartition, nSpokes, 1]);

    transformedCoordinates = reshape( ...
        reshape(coordinates4, [], 3) * coordinateTransform, ...
        [nReadout, nPartition, nSpokes, 3]);
    kz4 = transformedCoordinates(:, :, :, 3);
    kzReference = repmat(reshape(kz4(1, :, 1), [1 nPartition 1]), ...
        [nReadout, 1, nSpokes]);
    if any(abs(kz4(:) - kzReference(:)) > p.Results.KzTolerance)
        error('reconstructGARSOSRadialZSlice:NonseparableKz', ...
            'Transformed kz must be constant across readout and spoke for each partition.');
    end
    kz = reshape(kz4(1, :, 1), [nPartition, 1]);

    kernelSamplingDensity = p.Results.KernelSamplingDensity;
    kernelSamples = generateKaiserBesselKernel( ...
        imageSizeZ, gridSizeZ, kernelSizeZ, kernelSamplingDensity);
    kernelSamples2 = [fliplr(kernelSamples(2:end)), kernelSamples];
    kKernel = fftlattice(1 / kernelSamplingDensity, 0, numel(kernelSamples2));
    offsets = -floor(kernelSizeZ / 2):floor(kernelSizeZ / 2);

    kGridZ = kz .* gridSizeZ;
    iiZ = round(kGridZ);
    dkZ = iiZ - kGridZ;
    kernelWeights = interp1(kKernel, kernelSamples2, ...
        dkZ + reshape(offsets, 1, 1, []));
    if any(~isfinite(kernelWeights(:)))
        error('reconstructGARSOSRadialZSlice:KernelLookup', ...
            'The transformed kz coordinates fall outside the Kaiser-Bessel lookup range.');
    end

    grid4 = zeros(nReadout, gridSizeZ, nSpokes, nCoil, ...
        'like', sampledData);
    sampleLinear = reshape( ...
        (1:nReadout).' + (0:nSpokes-1) * nReadout, ...
        [nReadout, 1, nSpokes]);
    sampleLinear = repmat(sampleLinear, [1, nPartition, 1]);
    for iOffset = 1:numel(offsets)
        zIndexGrid = mod(iiZ + offsets(iOffset), gridSizeZ) + 1;
        gridLinear = sampleLinear + ...
            (zIndexGrid - 1) * nReadout;
        gridLinear = reshape(gridLinear, [], 1);
        weights = repmat(repelem( ...
            reshape(kernelWeights(:, :, iOffset), [], 1), nReadout), ...
            nSpokes, 1);
        weights = weights(:);
        for iCoil = 1:nCoil
            values = reshape(data4(:, :, :, iCoil), [], 1) .* weights;
            gridVector = accumarray(gridLinear, values, ...
                [nReadout * gridSizeZ * nSpokes, 1], @sum, 0);
            grid4(:, :, :, iCoil) = grid4(:, :, :, iCoil) + ...
                reshape(gridVector, [nReadout, gridSizeZ, nSpokes]);
        end
    end

    % GpuReconstructur20181008 applies ifftn without an FFT shift. Its
    % subsequent shiftAndResize call maps the edge-origin grid to a
    % centered image grid, which is reproduced below.
    imageGrid4 = ifft(grid4, [], 2);
    imageGrid2 = reshape(permute(imageGrid4, [2 1 3 4]), ...
        [gridSizeZ, nReadout * nSpokes * nCoil]);
    resized2 = shiftAndResize(imageGrid2, [], ...
        [imageSizeZ, nReadout * nSpokes * nCoil], 'ee', 'me');
    sliceData = reshape(resized2(zIndex, :), ...
        [nReadout, nSpokes, nCoil]);

    diagnostics = struct();
    diagnostics.nReadout = nReadout;
    diagnostics.nPartition = nPartition;
    diagnostics.nSpokes = nSpokes;
    diagnostics.nCoil = nCoil;
    diagnostics.transformedKz = kz;
    diagnostics.kGridZ = kGridZ;
    diagnostics.gridSizeZ = gridSizeZ;
    diagnostics.imageSizeZ = imageSizeZ;
    diagnostics.kernelSizeZ = kernelSizeZ;
    diagnostics.kernelSamplingDensity = kernelSamplingDensity;
    diagnostics.zIndex = zIndex;
    diagnostics.dataShape = '[readout, spoke, coil, z]';
    if p.Results.ReturnAllZ
        allZData = reshape(permute(reshape(resized2, ...
            [imageSizeZ, nReadout, nSpokes, nCoil]), [2 3 4 1]), ...
            [nReadout, nSpokes, nCoil, imageSizeZ]);
        diagnostics.allZData = allZData;
    end
end
