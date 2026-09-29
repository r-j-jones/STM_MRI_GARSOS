function fig = visualizeSTMEigenValues( eigenValues, plotDir )

if nargin<2, plotDir = []; end

fig = figure('Position',[570 82 989 407],'Color','w'); 
imagesc(utils.mdisp(abs(eigenValues)));
axis tight;
axis image;
colorbar; 
colormap gray;
clim([0 1])
title('Eigenvalues of G matrices (normalized)');
set(gca,'FontSize',15);

if ~isempty(plotDir)
    if ~exist(plotDir,'dir'), mkdir(plotDir); end
    outPlotPath = fullfile(plotDir,'G_matrix_eigenValues.png');
    print(fig,outPlotPath,'-dpng','-r300');
end

end