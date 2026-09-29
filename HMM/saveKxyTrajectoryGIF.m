function saveKxyTrajectoryGIF(kxy, framesToPlot, outputDirectory)
%SAVEKXYTRAJECTORYGIF Save selected trajectory frames as an animated GIF.
%
%   saveKxyTrajectoryGIF(kxy, framesToPlot, outputDirectory)
%
% Inputs
%   kxy             Complex array with dimensions:
%                   readout samples x spokes per frame x frames.
%                   real(kxy) is kx and imag(kxy) is ky.
%
%   framesToPlot    Vector specifying the frames and their animation order.
%
%   outputDirectory Directory in which the GIF will be saved.
%
% Output file
%   kxy_trajectories.gif
%
% Example
%   saveKxyTrajectoryGIF(kxy, 1:size(kxy,3), 'trajectory_gif');
%
% USAGE:
%  framesToPlot = 1:133;
%  outputDirectory = fullfile(pwd, 'trajectory_plots');
%  saveKxyTrajectoryGIF(kxy, framesToPlot, outputDirectory);

    narginchk(3, 3);

    validateattributes(kxy, {'numeric'}, ...
        {'nonempty', '3d'}, mfilename, 'kxy', 1);

    nFrames = size(kxy, 3);
    framesToPlot = validateFrames(framesToPlot, nFrames);

    if ~(ischar(outputDirectory) || ...
            (isstring(outputDirectory) && isscalar(outputDirectory)))
        error('saveKxyTrajectoryGIF:InvalidOutputDirectory', ...
            'outputDirectory must be a character vector or string scalar.');
    end

    outputDirectory = char(outputDirectory);

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    outputPath = fullfile(outputDirectory, 'kxy_trajectories.gif');
    delayTime = 0.15;  % Seconds per GIF frame

    [xLimits, yLimits] = getKxyLimits(kxy, framesToPlot);

    fig = figure( ...
        'Color', 'w', ...
        'Visible', 'off', ...
        'Position', [1265 295 496 413]);

    for iAnimationFrame = 1:numel(framesToPlot)
        iFrame = framesToPlot(iAnimationFrame);

        clf(fig);
        ax = axes(fig);
        hold(ax, 'on');

        for iSpoke = 1:size(kxy, 2)
            trajectory = kxy(:, iSpoke, iFrame);

            plot(ax, real(trajectory), imag(trajectory), ...
                '-', 'Color', 'k', 'LineWidth', 1.5);
        end

        xlabel(ax, 'k_x');
        ylabel(ax, 'k_y');
        title(ax, sprintf('Frame %d', iFrame));
        set(ax,'FontSize',15);

        xlim(ax, xLimits);
        ylim(ax, yLimits);
        axis(ax, 'equal');
        axis(ax, 'square');
        grid(ax, 'on');
        box(ax, 'on');

        drawnow;

        rgbImage = print(fig, '-RGBImage', '-r100');
        [indexedImage, colorMap] = rgb2ind(rgbImage, 256);

        if iAnimationFrame == 1
            imwrite(indexedImage, colorMap, outputPath, 'gif', ...
                'LoopCount', Inf, ...
                'DelayTime', delayTime);
        else
            imwrite(indexedImage, colorMap, outputPath, 'gif', ...
                'WriteMode', 'append', ...
                'DelayTime', delayTime);
        end
    end

    close(fig);
end


function frames = validateFrames(frames, nFrames)

    validateattributes(frames, {'numeric'}, ...
        {'vector', 'real', 'finite', 'integer', 'positive', 'nonempty'}, ...
        mfilename, 'framesToPlot', 2);

    frames = frames(:).';

    if any(frames > nFrames)
        error('saveKxyTrajectoryGIF:FrameOutOfRange', ...
            'framesToPlot must contain values between 1 and %d.', nFrames);
    end
end


function [xLimits, yLimits] = getKxyLimits(kxy, frames)

    selected = kxy(:, :, frames);
    kx = real(selected);
    ky = imag(selected);

    finitePoints = isfinite(kx) & isfinite(ky);

    if ~any(finitePoints(:))
        error('saveKxyTrajectoryGIF:NoFinitePoints', ...
            'The selected frames contain no finite trajectory points.');
    end

    kx = kx(finitePoints);
    ky = ky(finitePoints);

    xLimits = paddedLimits(min(kx), max(kx));
    yLimits = paddedLimits(min(ky), max(ky));
end


function limits = paddedLimits(minimumValue, maximumValue)

    rangeValue = maximumValue - minimumValue;

    if rangeValue == 0
        padding = max(1, abs(minimumValue)) * 0.05;
    else
        padding = 0.05 * rangeValue;
    end

    limits = [minimumValue - padding, maximumValue + padding];
end