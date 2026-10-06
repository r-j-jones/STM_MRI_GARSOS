function save_SpatialCoeff_images( stm_recon, outGifDir )

if ~exist(outGifDir,'dir'), mkdir(outGifDir); end

coefields = fieldnames(stm_recon);
coefields = coefields(contains(coefields,'coef'));

for cc=1:numel(coefields)
    curf = coefields{cc};
    if contains(curf,'tik')
        regtype='PCG-Tikh';
    else
        regtype='PCG';
    end
    curcoefs = stm_recon.(curf);
    coefarr = reshape(curcoefs, [stm_recon.imageSize stm_recon.nBasis]);
    coefarr = coefarr(:,:,end:-1:1);
    coefarr = abs(coefarr);
    origcoefarr = coefarr;

    for kk=1:size(coefarr,3)
        curdat = coefarr(:,:,kk);
        if strcmp(regtype,'PCG')
            scalval = prctile(curdat(:),66);
            curdat = curdat/scalval;
            curdat(curdat>1)=1;
        else
            scalval = prctile(curdat(:),99);
            curdat = curdat/scalval;
            curdat(curdat>1)=1;
        end
        coefarr(:,:,kk) = curdat;
    end

    fig = figure('color','w','position',[2367 121 678 577]); %[2359 121 686 734]);
    im(coefarr);
    ccll=clim;
    clim([0 ccll(2)]);
    colorbar;
    title( sprintf('%s Spatial Coeff Images\n(Indv Scaled)',regtype) ,'FontSize',20);
    
    outfilename = sprintf('%s_SpatCoef_images_montage_IndvScaled.png',regtype);
    print(fig, fullfile(outGifDir,outfilename), '-dpng','-r300');

    clf;
    im(origcoefarr);
    if strcmp(regtype,'PCG')
        scalval = prctile(origcoefarr(:,:,1),66,'all');
    else
        scalval = prctile(origcoefarr(:,:,1),99,'all');
    end
    clim([0 scalval]);
    colorbar;
    title( sprintf('%s Spatial Coeff Images\n(Global Scaled)',regtype) ,'FontSize',20);
    
    outfilename = sprintf('%s_SpatCoef_images_montage_GlobScaled.png',regtype);
    print(fig, fullfile(outGifDir,outfilename), '-dpng','-r300');

end

