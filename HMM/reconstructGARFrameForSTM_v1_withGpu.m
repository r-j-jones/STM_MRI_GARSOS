function [kCal, outputs, diagnostics] = reconstructGARFrameForSTM_v1( ...
    sampledData, sampleCoordinates, densityCompensation, spokeIndex, ...
    sensitivityMaps, imageSize, gridSize, nSpokesPerFrame, zSlice, ...
    calibrationSize, varargin)
%RECONSTRUCTGARFRAMEFORSTM Prepare 2-D STM calibration data from GAR-SOS.
%
% This wrapper preserves the existing 3-D GPU gridding implementation while
% preparing image-space slices and Cartesian calibration data for STM.
%
% Inputs:
%   sampledData          [Ns x Nc] complex corrected radial samples.
%   sampleCoordinates    [Ns x 3] normalized [kx ky kz] coordinates.
%   densityCompensation  [Ns x 1] radial density-compensation weights.
%   spokeIndex           [Ns x 1] zero- or one-based spoke identifier.
%   sensitivityMaps      [Nx x Ny x Nz x Nc] complex maps, or an object with
%                        an A property containing that array.
%   imageSize            [1 x 3] reconstructed volume dimensions.
%   gridSize             [1 x 3] oversampled Cartesian grid dimensions.
%   nSpokesPerFrame      positive integer number of consecutive spokes/frame.
%   zSlice               one-based image-space z index.
%   calibrationSize      [NxCal NyCal] central rectangular k-space size.
%
% Name/value options:
%   'KernelSize'         [1 x 3], default [7 7 7].
%   'OutputIsGpuArray'   logical, default false.
%   'CoordinateTransform' [3 x 3] HMM trajectory transform, default eye(3).
%   'TranslationVector'  [3 x 1] HMM translation vector for phase correction,
%                        default zeros(3,1).
%   'VisibleObjectSize'  [1 x 3] HMM visible object size, default imageSize.
%   'PhaseWidth'         scalar HMM phase-width normalization, default 1.
%   'SpatialSigma'       scalar HMM Gaussian weighting parameter, default 0.
%   'DensityScale'       additional scalar applied after HMM density scaling,
%                        default 1.
%   'ReturnDiagnostics'  logical, default false.
%
% Outputs:
%   kCal       [NxCal x NyCal x NFrames] complex Cartesian calibration data.
%   outputs    Structure containing internally coil-combined reconstructed
%              volumes, selected z slices, and Cartesian k-space.
%   diagnostics Optional structure with per-frame inputs and interface sizes.
%
% The existing GpuReconstructur20181008 class applies its validated
% gridding-kernel deapodization and coil-profile weighting. This wrapper
% applies the original HMM coordinate, phase, Gaussian, and density scaling
% once per frame, then combines all coils in one grid call.

    p = inputParser;
    addParameter(p, 'KernelSize', [7 7 7], @(x) isnumeric(x) && isvector(x) && numel(x) == 3);
    addParameter(p, 'OutputIsGpuArray', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'CoordinateTransform', eye(3), @(x) isnumeric(x) && isequal(size(x), [3 3]) && all(isfinite(x(:))));
    addParameter(p, 'TranslationVector', zeros(3, 1), @(x) isnumeric(x) && numel(x) == 3 && all(isfinite(x(:))));
    addParameter(p, 'VisibleObjectSize', [], @(x) isnumeric(x) && (isempty(x) || (isvector(x) && numel(x) == 3 && all(isfinite(x(:))))));
    addParameter(p, 'PhaseWidth', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x ~= 0);
    addParameter(p, 'SpatialSigma', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(p, 'DensityScale', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
    addParameter(p, 'ReturnDiagnostics', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    validateattributes(sampledData, {'numeric'}, {'2d', 'nonempty'}, mfilename, 'sampledData');
    validateattributes(sampleCoordinates, {'numeric'}, {'2d', 'ncols', 3}, mfilename, 'sampleCoordinates');
    validateattributes(densityCompensation, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'densityCompensation');
    validateattributes(spokeIndex, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'spokeIndex');
    validateattributes(imageSize, {'numeric'}, {'vector', 'numel', 3, 'integer', 'positive'}, mfilename, 'imageSize');
    validateattributes(gridSize, {'numeric'}, {'vector', 'numel', 3, 'integer', 'positive'}, mfilename, 'gridSize');
    validateattributes(nSpokesPerFrame, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'nSpokesPerFrame');
    validateattributes(zSlice, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'zSlice');
    validateattributes(calibrationSize, {'numeric'}, {'vector', 'numel', 2, 'integer', 'positive'}, mfilename, 'calibrationSize');

    [nSamples, nChannels] = size(sampledData);
    coordinateTransform = p.Results.CoordinateTransform;
    translationVector = p.Results.TranslationVector(:);
    visibleObjectSize = p.Results.VisibleObjectSize;
    if isempty(visibleObjectSize)
        visibleObjectSize = imageSize;
    end
    visibleObjectSize = visibleObjectSize(:).';
    if size(sampleCoordinates, 1) ~= nSamples
        error('reconstructGARFrameForSTM:SampleCountMismatch', ...
            'sampledData and sampleCoordinates must have the same number of rows.');
    end
    if numel(densityCompensation) ~= nSamples || numel(spokeIndex) ~= nSamples
        error('reconstructGARFrameForSTM:SampleCountMismatch', ...
            'densityCompensation and spokeIndex must contain one value per sample.');
    end

    sensitivityArray = localSensitivityArray(sensitivityMaps);
    if ndims(sensitivityArray) ~= 4
        error('reconstructGARFrameForSTM:SensitivityDimensions', ...
            'sensitivityMaps must contain [Nx Ny Nz Nc] data.');
    end
    if size(sensitivityArray, 4) ~= nChannels
        error('reconstructGARFrameForSTM:ChannelMismatch', ...
            'Sensitivity-map coil count must match sampledData.');
    end
    sensitivitySize = size(sensitivityArray);
    if any(sensitivitySize(1:3) ~= imageSize)
        error('reconstructGARFrameForSTM:SensitivityDimensions', ...
            'Sensitivity-map spatial dimensions must match imageSize.');
    end
    if zSlice > imageSize(3)
        error('reconstructGARFrameForSTM:InvalidSlice', ...
            'zSlice must be within the reconstructed image volume.');
    end
    if any(calibrationSize > imageSize(1:2))
        error('reconstructGARFrameForSTM:InvalidCalibrationSize', ...
            'calibrationSize cannot exceed the reconstructed in-plane size.');
    end
    if any(~isfinite(sampleCoordinates(:))) || any(~isfinite(densityCompensation(:))) || ...
            any(~isfinite(spokeIndex(:)))
        error('reconstructGARFrameForSTM:NonFiniteInput', ...
            'Coordinates, density compensation, and spoke indices must be finite.');
    end

    spokeValues = unique(spokeIndex(:), 'stable');
    nFrames = floor(numel(spokeValues) / nSpokesPerFrame);
    if nFrames < 2
        error('reconstructGARFrameForSTM:InsufficientSpokes', ...
            'At least two complete temporal frames are required for 3-D kCal.');
    end
    spokeValues = spokeValues(1:nFrames * nSpokesPerFrame);

    sensitivitySlice = sensitivityArray(:, :, zSlice, :);
    sensitivitySlice = reshape(sensitivitySlice, [imageSize(1:2), nChannels]);
    reconstructedVolumes = zeros([imageSize, nFrames], 'like', sampledData);
    senseImages = zeros([imageSize(1:2), nFrames], 'like', sampledData);
    kSpace2D = zeros([imageSize(1:2), nFrames], 'like', sampledData);
    gpuSensitivityArray = gpuArray(sensitivityArray);

    frameSampleIndices = cell(1, nFrames);
    frameSpokes = zeros(nSpokesPerFrame, nFrames);
    for iFrame = 1:nFrames
        frameSpokes(:, iFrame) = spokeValues((iFrame - 1) * nSpokesPerFrame + (1:nSpokesPerFrame));
        frameSampleIndices{iFrame} = ismember(spokeIndex, frameSpokes(:, iFrame));
    end

    dcf = densityCompensation(:) * p.Results.DensityScale;
    for iFrame = 1:nFrames
        sampleMask = frameSampleIndices{iFrame};
        frameData = sampledData(sampleMask, :);
        frameCoordinates = sampleCoordinates(sampleMask, :);
        frameDCF = dcf(sampleMask);

        if isempty(frameData)
            error('reconstructGARFrameForSTM:EmptyFrame', ...
                'Frame %d contains no samples.', iFrame);
        end

        frameTrajectory = frameCoordinates * coordinateTransform;
        phaseCorrection = exp(1i * 2 * pi * (frameCoordinates * translationVector));
        d2 = sum(frameTrajectory.^2, 2);
        gaussianWeight = exp(-2 * pi * d2 * p.Results.SpatialSigma.^2);
        densityWeight = frameDCF * prod(visibleObjectSize(1:2)) * 2 / p.Results.PhaseWidth;
        weightedFrameData = bsxfun(@times, frameData, ...
            densityWeight .* gaussianWeight .* phaseCorrection);

        gpuFrameTrajectory = gpuArray(frameTrajectory);
        gpuWeightedFrameData = gpuArray(weightedFrameData);
        reconstructor = GpuReconstructur20181008( ...
            gridSize, imageSize, p.Results.KernelSize, gpuFrameTrajectory, ...
            gpuSensitivityArray);
        reconstructedVolume = reconstructor.grid( ...
            gpuWeightedFrameData, p.Results.OutputIsGpuArray);

        if any(~isfinite(reconstructedVolume(:)))
            error('reconstructGARFrameForSTM:NonFiniteOutput', ...
                'Frame %d produced NaN or Inf values during gridding.', iFrame);
        end
        reconstructedVolumes(:, :, :, iFrame) = reconstructedVolume;
        senseImage = reconstructedVolume(:, :, zSlice);
        senseImages(:, :, iFrame) = senseImage;
        kSpace2D(:, :, iFrame) = fftshift(fft2(ifftshift(senseImage)));

        if any(~isfinite(senseImage(:)))
            error('reconstructGARFrameForSTM:NonFiniteOutput', ...
                'Frame %d produced NaN or Inf values during reconstruction.', iFrame);
        end

        if p.Results.ReturnDiagnostics
            diagnostics.frameData{iFrame} = frameData;
            diagnostics.frameTrajectory{iFrame} = frameTrajectory;
            diagnostics.frameDCF{iFrame} = frameDCF;
        end
    end

    kCal = localCentralRegion(kSpace2D, calibrationSize);
    if any(~isfinite(kCal(:)))
        error('reconstructGARFrameForSTM:NonFiniteOutput', ...
            'The extracted kCal contains NaN or Inf values.');
    end

    outputs = struct();
    outputs.reconstructedVolumes = reconstructedVolumes;
    outputs.sensitivitySlice = sensitivitySlice;
    outputs.senseImages = senseImages;
    outputs.kSpace2D = kSpace2D;

    if p.Results.ReturnDiagnostics
        diagnostics.nSamples = nSamples;
        diagnostics.nChannels = nChannels;
        diagnostics.nFrames = nFrames;
        diagnostics.frameSpokes = frameSpokes;
        diagnostics.frameSampleCounts = cellfun(@numel, frameSampleIndices);
        diagnostics.coordinateTransform = coordinateTransform;
        diagnostics.translationVector = translationVector;
        diagnostics.visibleObjectSize = visibleObjectSize;
        diagnostics.phaseWidth = p.Results.PhaseWidth;
        diagnostics.spatialSigma = p.Results.SpatialSigma;
        diagnostics.kCal = kCal;
    else
        diagnostics = struct();
    end
end

function sensitivityArray = localSensitivityArray(sensitivityMaps)
    if isstruct(sensitivityMaps) && isfield(sensitivityMaps, 'A')
        sensitivityArray = sensitivityMaps.A;
    elseif isobject(sensitivityMaps) && isprop(sensitivityMaps, 'A')
        sensitivityArray = sensitivityMaps.A;
    else
        sensitivityArray = sensitivityMaps;
    end
end

function kCal = localCentralRegion(kSpace2D, calibrationSize)
    imageSize = size(kSpace2D);
    center = floor(imageSize(1:2) / 2) + 1;
    startIndex = center - floor(calibrationSize / 2);
    index1 = startIndex(1) + (0:calibrationSize(1) - 1);
    index2 = startIndex(2) + (0:calibrationSize(2) - 1);
    kCal = kSpace2D(index1, index2, :);
end