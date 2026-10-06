function make_kSpace2D_gifs( kSpace2D, nSpokesPerFrame, displaySliceIndex, framesToUse, outGifDir )
%make_senseImages_gifs( senseImages, nSpokesPerFrame, displaySliceIndex, framesToUse )

if ~exist(outGifDir,'dir'), mkdir(outGifDir); end

% gifdata = abs(squeeze(senseImages(:,:,50,1:100)));
gifdata = abs(squeeze(kSpace2D(:,:,displaySliceIndex,framesToUse)));

int_im = prctile(abs(gifdata(:)),80);
gifbasename = sprintf('kSpace2D_%d-spokes_frames%d-%d_z-%d.gif', ...
    nSpokesPerFrame, framesToUse(1), framesToUse(end), displaySliceIndex);
gifname = fullfile(outGifDir,gifbasename);
im_g = utils.gif_image_garsos(gifname, gifdata, int_im, 0.15);

end

