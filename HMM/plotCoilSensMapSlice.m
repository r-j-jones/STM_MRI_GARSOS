function plotCoilSensMapSlice( senseMapsSlice, outputPlotDirectory )

coilFig = figure('position',[1984 272 699 568],'color','w');
h = im('row',4,'col',5,abs(senseMapsSlice));
title('Coil sensitivity maps slice');
set(gca,'FontSize',15);
% plasmamap = slanCM('plasma'); plasmamap = cat(1,[0 0 0],plasmamap);
% colormap(plasmamap);
colormap(gray);
colorbar;
print(coilFig,fullfile(outputPlotDirectory,'coil-sens-map-slice.png'),'-dpng','-r300');

end