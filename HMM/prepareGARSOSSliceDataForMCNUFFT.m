function prepared = prepareGARSOSSliceDataForMCNUFFT( ...
    sampledData, kSpaceLocations, spokeIndex, nReadout, nPartition, ...
    nSpokes, nbSpokesPerFrame, zIndex, varargin)
%prepareGARSOSSLICEDATAFORMCNUFFT Prepare slice-wise stack-of-stars data.
%
% Restores the flattened preparation output to:
%   [readout, partition, spoke, coil]
% and applies a centered unitary inverse FFT along the Cartesian partition
% dimension before selecting an image-space z slice.
%
% Inputs:
%   sampledData       [nReadout*nPartition*nSpokes, nCoil]
%   kSpaceLocations   [nReadout*nPartition*nSpokes, 3], [kx ky kz]
%   spokeIndex        one spoke identifier per flattened sample
%   nReadout          q.image.NCol
%   nPartition        number of partition samples in sampledData
%   nSpokes           number of spoke/line samples in sampledData
%   nbSpokesPerFrame  positive spokes per temporal frame
%   zIndex            one-based image-space z index
%
% Name/value:
%   'KzTransform'     'ifftshift-then-fftshift' (default) or
%                     'fftshift-then-ifftshift'
%   'CoordinateTransform' 3x3 trajectory transform, default eye(3)
%   'TranslationVector'   3x1 Fourier phase translation, default zeros(3,1)
%
% Outputs:
%   prepared.data     [nReadout, nbSpokesPerFrame, nCoil, nFrame]
%   prepared.kxy      [nReadout, nbSpokesPerFrame, nFrame] complex kx+i*ky
%   prepared.frameSpokes [nbSpokesPerFrame, nFrame]
%   prepared.nDiscardedSpokes
%   prepared.nDiscardedSamples
%   prepared.kSpaceLocations4D [nReadout, nPartition, nSpokes, 3]
%
% The input trajectory uses normalized cycles/FOV-like coordinates. The
% 2*pi conversion is intentionally deferred to the NUFFT-plan helper.

    p = inputParser;
    addParameter(p, 'KzTransform', 'ifftshift-then-fftshift', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'CoordinateTransform', eye(3), ...
        @(x) isnumeric(x) && isequal(size(x), [3 3]) && all(isfinite(x(:))));
    addParameter(p, 'TranslationVector', zeros(3, 1), ...
        @(x) isnumeric(x) && numel(x) == 3 && all(isfinite(x(:))));
    parse(p, varargin{:});
    kzTransform = char(p.Results.KzTransform);
    coordinateTransform = p.Results.CoordinateTransform;
    translationVector = p.Results.TranslationVector(:);

    validateattributes(sampledData, {'numeric'}, {'2d', 'nonempty'});
    validateattributes(kSpaceLocations, {'numeric'}, {'2d', 'ncols', 3});
    validateattributes(spokeIndex, {'numeric'}, {'vector', 'nonempty'});
    validateattributes(nReadout, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(nPartition, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(nSpokes, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(nbSpokesPerFrame, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(zIndex, {'numeric'}, {'scalar', 'integer', 'positive'});

    nRows = nReadout * nPartition * nSpokes;
    if size(sampledData, 1) ~= nRows || size(kSpaceLocations, 1) ~= nRows || ...
            numel(spokeIndex) ~= nRows
        error('prepareGARSOSSliceDataForMCNUFFT:SizeMismatch', ...
            'Flattened inputs must contain nReadout*nPartition*nSpokes rows.');
    end
    if zIndex > nPartition
        error('prepareGARSOSSliceDataForMCNUFFT:InvalidZIndex', ...
            'zIndex must be within 1:nPartition.');
    end
    if any(~isfinite(sampledData(:))) || any(~isfinite(kSpaceLocations(:))) || ...
            any(~isfinite(spokeIndex(:)))
        error('prepareGARSOSSliceDataForMCNUFFT:NonFiniteInput', ...
            'Data, trajectory, and spoke identifiers must be finite.');
    end

    nCoil = size(sampledData, 2);
    data4 = reshape(sampledData, [nReadout, nPartition, nSpokes, nCoil]);
    locations4 = reshape(kSpaceLocations, [nReadout, nPartition, nSpokes, 3]);
    spoke4 = reshape(spokeIndex, [nReadout, nPartition, nSpokes]);
    rawLocations4 = locations4;
    phase4 = exp(1i * 2 * pi * ...
        reshape(sum(rawLocations4 .* reshape(translationVector, [1 1 1 3]), 4), ...
        [nReadout, nPartition, nSpokes]));
    data4 = data4 .* reshape(phase4, [nReadout, nPartition, nSpokes, 1]);
    locations4 = reshape(reshape(rawLocations4, [], 3) * coordinateTransform, ...
        [nReadout, nPartition, nSpokes, 3]);

    spokeValues = reshape(spoke4(1, 1, :), [nSpokes, 1]);
    expectedSpokes = repmat(reshape(spokeValues, [1 1 nSpokes]), ...
        [nReadout, nPartition, 1]);
    if any(spoke4(:) ~= expectedSpokes(:))
        error('prepareGARSOSSliceDataForMCNUFFT:SpokeOrdering', ...
            'Each restored spoke block must have one constant spoke identifier.');
    end
    expectedKz = repmat(locations4(1, :, :, 3), [nReadout, 1, 1]);
    if ~isempty(find(locations4(:, :, :, 3) ~= expectedKz, 1))
        error('prepareGARSOSSliceDataForMCNUFFT:KzOrdering', ...
            'Each readout/spoke must have a consistent Cartesian kz coordinate.');
    end
    expectedKxy = repmat(locations4(1, 1, :, 1:2), [nReadout, nPartition, 1, 1]);
    if ~isempty(find(locations4(:, :, :, 1:2) ~= expectedKxy, 1))
        error('prepareGARSOSSliceDataForMCNUFFT:TrajectoryOrdering', ...
            'The in-plane trajectory must repeat across Cartesian partitions.');
    end

    nFrame = floor(nSpokes / nbSpokesPerFrame);
    nDiscardedSpokes = nSpokes - nFrame * nbSpokesPerFrame;
    nDiscardedSamples = nDiscardedSpokes * nReadout * nPartition;
    if nFrame < 1
        error('prepareGARSOSSliceDataForMCNUFFT:InsufficientSpokes', ...
            'No complete frame is available for the requested frame size.');
    end

    frameSpokes = reshape( ...
        spokeValues(1:nFrame * nbSpokesPerFrame), ...
        [nbSpokesPerFrame, nFrame]);

    switch kzTransform
        case 'ifftshift-then-fftshift'
            % Default: centered kz -> centered image-space z convention.
            dataZ4 = fftshift(ifft(ifftshift(data4, 2), [], 2), 2) * sqrt(nPartition);
            % Alternative convention to test if the z orientation is wrong:
            % dataZ4 = ifftshift(ifft(fftshift(data4, 2), [], 2), 2) * sqrt(nPartition);
        case 'fftshift-then-ifftshift'
            % Alternative: centered kz ordering with the opposite shift pair.
            dataZ4 = ifftshift(ifft(fftshift(data4, 2), [], 2), 2) * sqrt(nPartition);
        otherwise
            error('prepareGARSOSSliceDataForMCNUFFT:InvalidKzTransform', ...
                'Unknown KzTransform option: %s', kzTransform);
    end

    selectedZ = dataZ4(:, zIndex, 1:nFrame * nbSpokesPerFrame, :);
    selectedZ = reshape(selectedZ, ...
        [nReadout, nbSpokesPerFrame, nFrame, nCoil]);
    selectedZ = permute(selectedZ, [1 2 4 3]);

    kxy4 = locations4(:, :, :, 1) + 1i * locations4(:, :, :, 2);
    kxyFrame = kxy4(:, 1, 1:nFrame * nbSpokesPerFrame);
    kxyFrame = reshape(kxyFrame, [nReadout, nbSpokesPerFrame, nFrame]);

    prepared = struct();
    prepared.data = selectedZ;
    prepared.kxy = kxyFrame;
    prepared.frameSpokes = frameSpokes;
    prepared.nDiscardedSpokes = nDiscardedSpokes;
    prepared.nDiscardedSamples = nDiscardedSamples;
    prepared.nReadout = nReadout;
    prepared.nPartition = nPartition;
    prepared.nSpokes = nSpokes;
    prepared.nbSpokesPerFrame = nbSpokesPerFrame;
    prepared.nCoil = nCoil;
    prepared.nFrame = nFrame;
    prepared.zIndex = zIndex;
    prepared.kzTransform = kzTransform;
    prepared.coordinateTransform = coordinateTransform;
    prepared.translationVector = translationVector;
    prepared.rawKSpaceLocations4D = rawLocations4;
    prepared.phaseCorrection = phase4;
    prepared.dataShape = '[readout, spoke, coil, frame]';
    prepared.kxyShape = '[readout, spoke, frame]';
    prepared.kSpaceLocations4D = locations4;
end
