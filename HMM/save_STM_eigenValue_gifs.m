function save_STM_eigenValue_gifs( result, nSpokesPerFrame, outGifDir, int_im, delayTime )

if nargin<4 || isempty(int_im), int_im = 0.2; end
if nargin<5 || isempty(delayTime), delayTime = 0.15; end

if ~exist(outGifDir,'dir'), mkdir(outGifDir); end


nStms = size(result.ST_maps,4);
for jj=1:nStms
    gifdata = abs(squeeze(result.ST_maps(:,:,:,jj)));
    ffname = sprintf('st_map_%d-spokes_%d.gif',nSpokesPerFrame,jj);
    gifname = fullfile(outGifDir,ffname); %['/RadOnc-MRI1/Student_Folder/rjones/STM/InitialTest_09-20-2026/' ffname];
    im_g = utils.gif_image_garsos(gifname, gifdata, int_im, delayTime);
end

gifdata = abs(squeeze(result.eigenValues(:,:,:)));
ffname = sprintf('eigenvals_%d-spokes.gif',nSpokesPerFrame);
gifname = fullfile(outGifDir,ffname); %['/RadOnc-MRI1/Student_Folder/rjones/STM/InitialTest_09-20-2026/' ffname];
im_g = utils.gif_image_garsos(gifname, gifdata, int_im, delayTime);

end
