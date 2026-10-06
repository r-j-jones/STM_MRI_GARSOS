function prepared = prepareGARSOSSliceDataForMCNUFFT( ...
    sampledData, kSpaceLocations, spokeIndex, nReadout, nPartition, ...
    nSpokes, nbSpokesPerFrame, zIndex, varargin)
%prepareGARSOSSLICEDATAFORMCNUFFT Prepare slice-wise stack-of-stars data.
%
% Restores the flattened preparation output to:
%   [readout, partition, spoke, coil]
% and applies either the legacy partition IFFT or the full-reconstruction
% compatible gridded z-only reconstruction before selecting a z slice.
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
%   'ZReconstructionMethod'  'direct-ifft' (default) or 'gridded'
%   'ImageSizeZ'         target image-space z size for 'gridded'
%   'GridSizeZ'          oversampled Cartesian z size for 'gridded'
%   'KernelSizeZ'        Kaiser-Bessel z kernel width for 'gridded'
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
    addParameter(p, 'ZReconstructionMethod', 'direct-ifft', ...
        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'ImageSizeZ', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0 && mod(x, 1) == 0));
    addParameter(p, 'GridSizeZ', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0 && mod(x, 1) == 0));
    addParameter(p, 'KernelSizeZ', 7, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0 && mod(x, 2) == 1);
    parse(p, varargin{:});
    kzTransform = char(p.Results.KzTransform);
    coordinateTransform = p.Results.CoordinateTransform;
    translationVector = p.Results.TranslationVector(:);
    zReconstructionMethod = lower(char(p.Results.ZReconstructionMethod));

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
    if strcmp(zReconstructionMethod, 'direct-ifft') && zIndex > nPartition
        error('prepareGARSOSSliceDataForMCNUFFT:InvalidZIndex', ...
            'zIndex must be within 1:nPartition.');
    end
    if ~ismember(zReconstructionMethod, {'direct-ifft', 'gridded'})
        error('prepareGARSOSSliceDataForMCNUFFT:InvalidZReconstructionMethod', ...
            'ZReconstructionMethod must be ''direct-ifft'' or ''gridded''.');
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
    expectedKxy = repmat(locations4(:, 1, :, 1:2), [1, nPartition, 1, 1]);
    if any(locations4(:, :, :, 1:2) ~= expectedKxy, 'all')
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

    switch zReconstructionMethod
        case 'direct-ifft'
            switch kzTransform
                case 'ifftshift-then-fftshift'
                    % Legacy centered partition-FFT convention.
                    dataZ4 = fftshift(ifft(ifftshift(data4, 2), [], 2), 2) * sqrt(nPartition);
                case 'fftshift-then-ifftshift'
                    dataZ4 = ifftshift(ifft(fftshift(data4, 2), [], 2), 2) * sqrt(nPartition);
                otherwise
                    error('prepareGARSOSSliceDataForMCNUFFT:InvalidKzTransform', ...
                        'Unknown KzTransform option: %s', kzTransform);
            end
            selectedZ = dataZ4(:, zIndex, 1:nFrame * nbSpokesPerFrame, :);
        case 'gridded'
            imageSizeZ = p.Results.ImageSizeZ;
            gridSizeZ = p.Results.GridSizeZ;
            if isempty(imageSizeZ) || isempty(gridSizeZ)
                error('prepareGARSOSSliceDataForMCNUFFT:MissingZGridSize', ...
                    'ImageSizeZ and GridSizeZ are required for gridded reconstruction.');
            end
            [selectedZ, zDiagnostics] = reconstructGARSOSRadialZSlice( ...
                sampledData, kSpaceLocations, nReadout, nPartition, nSpokes, ...
                coordinateTransform, imageSizeZ, gridSizeZ, ...
                p.Results.KernelSizeZ, zIndex, ...
                'TranslationVector', translationVector, ...
                'ReturnAllZ', false);
            selectedZ = reshape(selectedZ(:, 1:nFrame * nbSpokesPerFrame, :), ...
                [nReadout, nFrame * nbSpokesPerFrame, nCoil]);
            selectedZ = reshape(selectedZ, ...
                [nReadout, nbSpokesPerFrame, nFrame, nCoil]);
            selectedZ = permute(selectedZ, [1 2 4 3]);
    end

    kxy4 = locations4(:, :, :, 1) + 1i * locations4(:, :, :, 2);
    kxyFrame = kxy4(:, 1, 1:nFrame * nbSpokesPerFrame);
    kxyFrame = reshape(kxyFrame, [nReadout, nbSpokesPerFrame, nFrame]);

    prepared = makePreparedOutput(selectedZ, kxyFrame, frameSpokes, ...
        'nDiscardedSpokes', nDiscardedSpokes, ...
        'nDiscardedSamples', nDiscardedSamples, ...
        'nReadout', nReadout, ...
        'nPartition', nPartition, ...
        'nSpokes', nSpokes, ...
        'nbSpokesPerFrame', nbSpokesPerFrame, ...
        'nCoil', nCoil, ...
        'nFrame', nFrame, ...
        'zIndex', zIndex, ...
        'kzTransform', kzTransform, ...
        'zReconstructionMethod', zReconstructionMethod, ...
        'coordinateTransform', coordinateTransform, ...
        'translationVector', translationVector, ...
        'rawKSpaceLocations4D', rawLocations4, ...
        'phaseCorrection', phase4, ...
        'dataShape', '[readout, spoke, coil, frame]', ...
        'kxyShape', '[readout, spoke, frame]', ...
        'kSpaceLocations4D', locations4);
    if strcmp(zReconstructionMethod, 'gridded')
        prepared.zDiagnostics = zDiagnostics;
    end
end

function prepared = makePreparedOutput(data, kxy, frameSpokes, varargin)
%makePreparedOutput Assemble the prepared-data output structure.
    if mod(numel(varargin), 2) ~= 0
        error('prepareGARSOSSliceDataForMCNUFFT:InvalidOutputOptions', ...
            'Output options must be supplied as name/value pairs.');
    end

    prepared = struct('data', data, 'kxy', kxy, 'frameSpokes', frameSpokes);
    for optionIndex = 1:2:numel(varargin)
        optionName = varargin{optionIndex};
        if ~(ischar(optionName) || ...
                (isstring(optionName) && isscalar(optionName)))
            error('prepareGARSOSSliceDataForMCNUFFT:InvalidOutputOptionName', ...
                'Output option names must be text.');
        end
        prepared.(char(optionName)) = varargin{optionIndex + 1};
    end
end
