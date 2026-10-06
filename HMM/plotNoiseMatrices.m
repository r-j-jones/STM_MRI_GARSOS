function figs = plotNoiseMatrices(debugStructure)
%PLOTNOISEMATRICES Create separate real/imaginary/magnitude figures.

if exist('viridis', 'file') ~= 2
    addpath('/path/to/viridis');
end

figs = [];
figs.noisecovariance = plotComplexMatrix(debugStructure.noiseCovarianceMatrix, 'Noise covariance matrix');
figs.noisewhitening = plotComplexMatrix(debugStructure.noiseWhiteningTransform, 'Noise whitening transform');

end

function fig = plotComplexMatrix(A, figureName)

fig = figure('Name', figureName, 'Color','w','Position',[692 505 1113 352]);

subplot(1, 3, 1);
imagesc(real(A));
axis image;
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Real');
set(gca,'FontSize',14);

subplot(1, 3, 2);
imagesc(imag(A));
axis image;
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Imaginary');
set(gca,'FontSize',14);

subplot(1, 3, 3);
imagesc(abs(A));
axis image;
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Magnitude');
set(gca,'FontSize',14);

sgtitle(figureName,'FontSize',18);
colormap gray;
colormap(viridis);

end