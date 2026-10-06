function im_g = gif_image_garsos(filename_GIM, data_temp, int_im, delay_time)

Nt = size(data_temp, 3);

for idx = 1:Nt
   im_g = abs(data_temp(:,:,idx));
   im_g = abs(im_g)/int_im;
   im_rgb = im_g(:,:,[1 1 1]);
   
   im_rgb = insertText(im_rgb,[-5 -5],int2str(idx),...
       'BoxColor','black', ...
       'TextColor', 'r', ...
       'FontSize',40, ...
       'BoxOpacity', 0, ...
       'AnchorPoint', 'LeftTop');
   % 'Font', "Arial Unicode", ...
       

   [Ag,map] = rgb2ind(im_rgb,256);
    if idx == 1
        imwrite(Ag,map,filename_GIM,"gif","LoopCount",Inf,"DelayTime",delay_time);
    else
        imwrite(Ag,map,filename_GIM,"gif","WriteMode","append","DelayTime",delay_time);
    end

%     if idx == M
%         imwrite(im_g, ['POGM_tp_iter_' int2str(M) '.png']);
%     end

end