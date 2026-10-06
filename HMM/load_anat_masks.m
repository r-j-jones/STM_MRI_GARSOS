function arrvol = load_anat_masks( inputPath )

arrvol=[];
bodymaskfile = fullfile(inputPath,'anat-masks','body-mask.nii');
if exist(bodymaskfile,'file')
    arrvol = loadArrayVolumeFromNii(bodymaskfile);
end

