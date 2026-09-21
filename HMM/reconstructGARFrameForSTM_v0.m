function [kCal, outputs, diagnostics] = reconstructGARFrameForSTM_v0( ...
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
%   'DensityScale'       scalar applied after densityCompensation, default 1.
%   'SensitivityTolerance' scalar denominator tolerance, default 0.
%   'ReturnDiagnostics'  logical, default false.
%
% Outputs:
%   kCal       [NxCal x NyCal x NFrames] complex Cartesian calibration data.
%   outputs    Structure containing reconstructed coil volumes, z slices,
%              SENSE images, and Cartesian k-space.
%   diagnostics Optional structure with per-frame inputs and interface sizes.
%
% The existing GpuReconstructur20181008 class applies its validated
% gridding-kernel deapodization and coil-profile weighting. This wrapper
% does not replace or disable that behavior. Since each frame has its own
% trajectory subset and the class combines channels in one grid call, one
% reconstructor is created per coil and frame.

    p = inputParser;
    addParameter(p, 'KernelSize', [7 7 7], @(x) isnumeric(x) && isvector(x) && numel(x) == 3);
    addParameter(p, 'OutputIsGpuArray', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'DensityScale', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
    addParameter(p, 'SensitivityTolerance', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
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
    denominator = sum(abs(sensitivitySlice).^2, 3);
    if p.Results.SensitivityTolerance > 0
        denominator(denominator < p.Results.SensitivityTolerance) = 0;
    end

    coilVolumes = zeros([imageSize, nChannels, nFrames], 'like', sampledData);
    zSlices = zeros([imageSize(1:2), nChannels, nFrames], 'like', sampledData);
    senseImages = zeros([imageSize(1:2), nFrames], 'like', sampledData);
    kSpace2D = zeros([imageSize(1:2), nFrames], 'like', sampledData);

    frameSampleIndices = cell(1, nFrames);
    frameSpokes = zeros(nSpokesPerFrame, nFrames);
    for iFrame = 1:nFrames
        frameSpokes(:, iFrame) = spokeValues((iFrame - 1) * nSpokesPerFrame + (1:nSpokesPerFrame));
        frameSampleIndices{iFrame} = ismember(spokeIndex, frameSpokes(:, iFrame));
    end

    dcf = densityCompensation(:) * p.Results.DensityScale;
    for iFrame = 101:105 %1:nFrames
        tic
        sampleMask = frameSampleIndices{iFrame};
        frameData = sampledData(sampleMask, :);
        frameTrajectory = sampleCoordinates(sampleMask, :);
        frameDCF = dcf(sampleMask);

        if isempty(frameData)
            error('reconstructGARFrameForSTM:EmptyFrame', ...
                'Frame %d contains no samples.', iFrame);
        end

        weightedFrameData = bsxfun(@times, frameData, frameDCF);

        %[ OPTION 1:
        %[ Pass all coils to reconstructor
        reconstructor = GpuReconstructur20181008( ...
            gridSize, imageSize, p.Results.KernelSize, frameTrajectory, ...
            sensitivityArray(:, :, :, :));
        %[ Recon gridded coil combined volume
        senseImageVol = reconstructor.grid( ...
            weightedFrameData(:, :), p.Results.OutputIsGpuArray);
        senseImage = senseImageVol(:,:,zSlice);

        % %[ OPTION 2:
        % %[ Pass each coil separately to reconstructor
        % for iChannel = 1:nChannels
        %     reconstructor = GpuReconstructur20181008( ...
        %         gridSize, imageSize, p.Results.KernelSize, frameTrajectory, ...
        %         sensitivityArray(:, :, :, iChannel));
        %     coilVolume = reconstructor.grid( ...
        %         weightedFrameData(:, iChannel), p.Results.OutputIsGpuArray);
        %     if any(~isfinite(coilVolume(:)))
        %         error('reconstructGARFrameForSTM:NonFiniteOutput', ...
        %             'Frame %d, coil %d produced NaN or Inf values.', iFrame, iChannel);
        %     end
        %     coilVolumes(:, :, :, iChannel, iFrame) = coilVolume;
        %     zSlices(:, :, iChannel, iFrame) = coilVolume(:, :, zSlice);
        % end
        % %[ SENSE coil combination
        % denominatorForCombination = denominator;
        % denominatorForCombination(denominatorForCombination == 0) = 1;
        % senseImage = sum(conj(sensitivitySlice) .* zSlices(:, :, :, iFrame), 3) ./ ...
        %     denominatorForCombination;
        % senseImage(denominator == 0) = 0;

        % Store the reconstructed senseImage and the 2d kspace for zSlice
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
        toc
    end

    kCal = localCentralRegion(kSpace2D, calibrationSize);
    if any(~isfinite(kCal(:)))
        error('reconstructGARFrameForSTM:NonFiniteOutput', ...
            'The extracted kCal contains NaN or Inf values.');
    end

    outputs = struct();
    outputs.coilVolumes = coilVolumes;
    outputs.zSlices = zSlices;
    outputs.sensitivitySlice = sensitivitySlice;
    outputs.senseImages = senseImages;
    outputs.kSpace2D = kSpace2D;

    if p.Results.ReturnDiagnostics
        diagnostics.nSamples = nSamples;
        diagnostics.nChannels = nChannels;
        diagnostics.nFrames = nFrames;
        diagnostics.frameSpokes = frameSpokes;
        diagnostics.frameSampleCounts = cellfun(@numel, frameSampleIndices);
        diagnostics.senseDenominator = denominator;
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
