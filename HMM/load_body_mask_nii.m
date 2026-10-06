function arrvol = load_body_mask_nii( inputPath )

arrvol=[];
bodymaskfile = fullfile(inputPath,'anat-masks','body-mask.nii');
if exist(bodymaskfile,'file')
    arrvol = loadArrayVolumeFromNii(bodymaskfile);
end

