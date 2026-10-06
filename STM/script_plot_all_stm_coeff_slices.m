
% Plot STM spatial coeff maps

addpath('/RadOnc-MRI1/Student_Folder/rjones/STM/STM_MRI_GARSOS/STM');

plotdir = '/RadOnc-MRI1/Student_Folder/rjones/STM/InitialTest_09-20-2026/Comparison-Plots';

for ii=1:L
    ii_ = L-ii+1;
    fig = fimshow( abs(zmat(:,:,ii_)), [], true, sprintf('Coeff %d',ii) );

    outfile = fullfile(plotdir,sprintf('demo-cart2d_spat-coeff-%d.png',ii_));
    print(fig, outfile,'-dpng','-r300');

end




plotdir = '/RadOnc-MRI1/Student_Folder/rjones/STM/InitialTest_09-20-2026/Comparison-Plots';

for ii=1:L
    ii_ = L-ii+1;
    fig = fimshow( abs(zmat(:,:,ii_)), [], true, sprintf('Coeff %d',ii) );

    outfile = fullfile(plotdir,sprintf('demo-cart2d_spat-coeff-%d.png',ii_));
    print(fig, outfile,'-dpng','-r300');

end



