function visualize_kdata_cal_movie( kCal, zIndex )

imageSize2D = [224 224];

if nargin < 2
    zIndex = 50;      % selected z slice
end

figure;

for frameIndex = 1:size(kCal, 4)
    kSmall = kCal(:, :, zIndex, frameIndex);

    kPadded = zeros(imageSize2D, 'like', kSmall);

    rowStart = floor(imageSize2D(1) / 2) - floor(size(kSmall, 1) / 2) + 1;
    colStart = floor(imageSize2D(2) / 2) - floor(size(kSmall, 2) / 2) + 1;

    rows = rowStart:(rowStart + size(kSmall, 1) - 1);
    cols = colStart:(colStart + size(kSmall, 2) - 1);

    kPadded(rows, cols) = kSmall;

    imageInterpolated = fftshift(ifft2(ifftshift(kPadded)));
    imageInterpolated = imageInterpolated * ...
        prod(imageSize2D) / prod(size(kSmall));

    imagesc(abs(imageInterpolated));
    axis image off;
    colormap gray;
    colorbar;
    title(sprintf('Slice %d, frame %d', zIndex, frameIndex));
    drawnow;
    pause(0.4);
end