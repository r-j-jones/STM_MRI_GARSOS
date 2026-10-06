function outputs = make_stm_recon_intensity_plots( curdata, outGifDir, reconType, sliceNb, frameNb )

    outputs = [];

    nFrames = size(curdata,3);
    if nargin<5 || isempty(frameNb), frameNb=67; end
    if nargin<4 || isempty(sliceNb), sliceNb=-99; end

    % Crop to center roi rectangle
    centerRoiX = 41:168;   % Row indices
    centerRoiY = 50:175;   % Column indices

    % Reshape to [number of ROI pixels, number of slices]
    stmReconCrop = abs(curdata(centerRoiX, centerRoiY, :));
    roiBySlice = reshape(stmReconCrop, [], nFrames);


    %% Calculate roi stats
    roiMin    = min(roiBySlice, [], 1, 'omitnan');
    roiMean   = mean(roiBySlice, 1, 'omitnan');
    roiMedian = median(roiBySlice, 1, 'omitnan');
    roiMax    = max(roiBySlice, [], 1, 'omitnan');
    roiStd    = std(roiBySlice, 0, 1, 'omitnan');
    roiP05    = prctile(roiBySlice, 5, 1);   % NaN handling depends on MATLAB version
    roiP95    = prctile(roiBySlice, 95, 1);
    

    %% PCG
    
    % curdata = stm_recon.sense_recon_stm;

    %% Display ROI on selected slice
   
    roifig = figure('position',[1065 521 435 350],'color','w');
    imagesc(abs(curdata(:, :, frameNb)));
    axis image;
    colormap gray;
    colorbar;
    hold on;
    
    % Overlay crop ROI
    rectangle( ...
        'Position', [centerRoiY(1)-0.5, centerRoiX(1)-0.5, ...
                     numel(centerRoiY), numel(centerRoiX)], ...
        'EdgeColor', 'r', ...
        'LineWidth', 2);
    
    title(sprintf('%s Recon, frame %d', upper(reconType),frameNb));
    hold off;
    set(gca,'FontSize',14);
    ccll=clim;
    if contains(lower(reconType),'tik')
        clim([ccll(1) mean(roiP95)]);
    else
        clim([ccll(1) roiMean(frameNb)*5]);
    end
    drawnow; shg;
    outRoiFig = fullfile(outGifDir,sprintf('%s_roi-location_z-%d.png',reconType,sliceNb));
    print(roifig,outRoiFig,'-dpng','-r300');
    outputs.roiLocationFile = outRoiFig;


    %% Plot ROI stats

    



    frameNumbers = 1:nFrames;

    %% Plot 1: percentiles and range
    lfig2 = figure('color','w','position',[592 606 921 239]);
    plot(frameNumbers, roiMin,    'LineWidth', 1.5); hold on;
    plot(frameNumbers, roiP05,    'LineWidth', 1.5);
    plot(frameNumbers, roiMedian, 'LineWidth', 2);
    plot(frameNumbers, roiP95,    'LineWidth', 1.5);
    plot(frameNumbers, roiMax,    'LineWidth', 1.5);
    hold off;
    
    grid on;
    xlabel('Frame');
    ylabel('ROI intensity (magn)');
    title(sprintf('%s ROI Intensity Stats',upper(reconType)));
    legend('Min', '5th prctile', 'Median', ...
           '95th prctile', 'Max', ...
           'Location', 'eastoutside');
    set(gca,'FontSize',13);
    xlim([1 nFrames]);

    drawnow; shg;
    outLineplotFig2 = fullfile(outGifDir,sprintf('%s_roi-frame-stats_2_z-%d.png',reconType,sliceNb));
    print(lfig2,outLineplotFig2,'-dpng','-r300');
    outputs.roiStatsLine2File = outLineplotFig2;

    yyll = ylim;

    ylim([yyll(1) max(roiP95)*1.01]);
    shg;
    outLineplotFig2 = fullfile(outGifDir,sprintf('%s_roi-frame-stats_2_no-max_z-%d.png',reconType,sliceNb));
    print(lfig2,outLineplotFig2,'-dpng','-r300');
    outputs.roiStatsLine2File_nomax = outLineplotFig2;

    ylim([yyll(1) max(roiMedian)*1.01]);
    shg;
    outLineplotFig2 = fullfile(outGifDir,sprintf('%s_roi-frame-stats_2_no-max-p95_z-%d.png',reconType,sliceNb));
    print(lfig2,outLineplotFig2,'-dpng','-r300');
    outputs.roiStatsLine2File_nomaxp95 = outLineplotFig2;

    
    %% Plot 2: mean with ±1 standard-deviation shaded region
    lowerBound = roiMean - roiStd;
    upperBound = roiMean + roiStd;
    
    lfig1 = figure('color','w','position',[592 606 921 239]);
    fill([frameNumbers, fliplr(frameNumbers)], ...
         [lowerBound, fliplr(upperBound)], ...
         0.7*[0.75 0.85 1.00], ...
         'EdgeColor', 'none', ...
         'FaceAlpha', 0.4);
    hold on;
    
    plot(frameNumbers, roiMean, 'b-', 'LineWidth', 2);
    
    % hold off;
    % grid on;
    % xlabel('Frame');
    % ylabel('ROI magnitude');
    % title(sprintf('%s Mean \pm 1 Standard Deviation',reconType));
    % legend('\pm1 standard deviation', 'Mean', 'Location', 'best');

    grid on;
    xlabel('Frame');
    ylabel('ROI intensity (magn)');
    title(sprintf('%s ROI Mean +/- SD',upper(reconType)));
    legend('\pm1 SD', 'Mean',  ...
           'Location', 'eastoutside');
    set(gca,'FontSize',13);
    xlim([1 nFrames]);

    drawnow; shg;
    outLineplotFig1 = fullfile(outGifDir,sprintf('%s_roi-frame-stats_1_z-%d.png',reconType,sliceNb));
    print(lfig1,outLineplotFig1,'-dpng','-r300');
    outputs.roiStatsLine1File = outLineplotFig1;



    %% Make gif of stm recons

    delay_time = 0.15;
    int_im = mean(roiP95);

    % PCG
    gifFile = fullfile(outGifDir,sprintf('%s_250-iters_z-%d.gif',reconType,sliceNb));
    gifData = rot90( abs(curdata), -1 );

    im_g = utils.gif_image_garsos(gifFile, gifData, int_im, delay_time);
    outputs.stmReconImagesGif = gifFile;

end