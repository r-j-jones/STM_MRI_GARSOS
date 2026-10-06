function make_kCal_gifs( kCal, nSpokesPerFrame, displaySliceIndex, framesToUse, calibrationSize, outGifDir )
%make_senseImages_gifs( senseImages, nSpokesPerFrame, displaySliceIndex, framesToUse )

if ~exist(outGifDir,'dir'), mkdir(outGifDir); end

% gifdata = abs(squeeze(senseImages(:,:,50,1:100)));
gifdata = abs(squeeze(kCal(:,:,displaySliceIndex,framesToUse)));

int_im = prctile(abs(gifdata(:)),80);
gifbasename = sprintf('kCal_len-%d_%d-spokes_frames%d-%d_z-%d.gif', ...
    calibrationSize, nSpokesPerFrame, framesToUse(1), framesToUse(end), displaySliceIndex);
gifname = fullfile(outGifDir,gifbasename);
im_g = utils.gif_image_garsos(gifname, gifdata, int_im, 0.15);

end

