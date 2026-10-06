function make_senseImages_gifs( senseImages, nSpokesPerFrame, displaySliceIndex, framesToUse, outGifDir )
%make_senseImages_gifs( senseImages, nSpokesPerFrame, displaySliceIndex, framesToUse )

if ~exist(outGifDir,'dir'), mkdir(outGifDir); end

% gifdata = abs(squeeze(senseImages(:,:,50,1:100)));
gifdata = abs(squeeze(senseImages(:,:,displaySliceIndex,framesToUse)));

gifbasename = sprintf('senseImages_%d-spokes_frames%d-%d_z-%d.gif', ...
    nSpokesPerFrame, framesToUse(1), framesToUse(end), displaySliceIndex);
gifname = fullfile(outGifDir,gifbasename);
im_g = utils.gif_image_garsos(gifname, gifdata, 5, 0.15);

end

