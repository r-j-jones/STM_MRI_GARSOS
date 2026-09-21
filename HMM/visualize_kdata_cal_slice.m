function visualize_kdata_cal_slice( kCal, zIndex, frameIndex )

imageSize2D = [224 224];

if nargin < 2
    zIndex = 50;      % selected z slice
end
if nargin < 3
    frameIndex = 50;   % selected temporal frame
end

kSmall = kCal(:, :, zIndex, frameIndex);

kPadded = zeros(imageSize2D, 'like', kSmall);

smallCenter = floor(size(kSmall) / 2) + 1;
largeCenter = floor(imageSize2D / 2) + 1;

rowStart = largeCenter(1) - floor(size(kSmall, 1) / 2);
colStart = largeCenter(2) - floor(size(kSmall, 2) / 2);

rows = rowStart:(rowStart + size(kSmall, 1) - 1);
cols = colStart:(colStart + size(kSmall, 2) - 1);

kPadded(rows, cols) = kSmall;

imageInterpolated = fftshift(ifft2(ifftshift(kPadded)));

% Compensate for MATLAB's larger inverse-FFT normalization
imageInterpolated = imageInterpolated * ...
    prod(imageSize2D) / prod(size(kSmall));

figure;
imagesc(abs(imageInterpolated));
axis image off;
colormap gray;
colorbar;
title(sprintf('Slice %d, frame %d', zIndex, frameIndex));