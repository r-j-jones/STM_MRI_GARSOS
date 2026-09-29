function saveKxyFramePNGs(kxy, framesToPlot, outputDirectory)
%SAVEKXYFRAMEPNGS Plot selected radial trajectories and save them as PNGs.
%
%   saveKxyFramePNGs(kxy, framesToPlot, outputDirectory)
%
% Inputs
%   kxy             Complex array with dimensions:
%                   readout samples x spokes per frame x frames.
%                   real(kxy) is kx and imag(kxy) is ky.
%
%   framesToPlot    Vector containing the frame numbers to plot.
%
%   outputDirectory Directory in which PNG files will be saved.
%
% Example
%   saveKxyFramePNGs(kxy, [1 10 20], 'trajectory_pngs');
%
% USAGE:
%  framesToPlot = 1:133;
%  outputDirectory = fullfile(pwd, 'trajectory_plots');
%  saveKxyFramePNGs(kxy, framesToPlot, outputDirectory);

    narginchk(3, 3);

    validateattributes(kxy, {'numeric'}, ...
        {'nonempty', '3d'}, mfilename, 'kxy', 1);

    nFrames = size(kxy, 3);
    framesToPlot = validateFrames(framesToPlot, nFrames);

    if ~(ischar(outputDirectory) || ...
            (isstring(outputDirectory) && isscalar(outputDirectory)))
        error('saveKxyFramePNGs:InvalidOutputDirectory', ...
            'outputDirectory must be a character vector or string scalar.');
    end

    outputDirectory = char(outputDirectory);

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    [xLimits, yLimits] = getKxyLimits(kxy, framesToPlot);

    fig = figure('Color', 'w', 'Visible', 'off');

    for iFrame = framesToPlot
        clf(fig);
        ax = axes(fig);
        hold(ax, 'on');

        for iSpoke = 1:size(kxy, 2)
            trajectory = kxy(:, iSpoke, iFrame);

            plot(ax, real(trajectory), imag(trajectory), ...
                '-o', 'LineWidth', 1, 'Color','k', 'MarkerSize', 4);
        end

        xlabel(ax, 'k_x');
        ylabel(ax, 'k_y');
        title(ax, sprintf('Frame %d', iFrame));
        set(ax,'FontSize',15);

        xlim(ax, xLimits);
        ylim(ax, yLimits);
        axis(ax, 'equal');
        grid(ax, 'on');
        box(ax, 'on');

        outputPath = fullfile(outputDirectory, ...
            sprintf('kxy_frame_%04d.png', iFrame));

        figure(fig); drawnow;
        print(gcf, outputPath, '-dpng', '-r300');
    end

    close(fig);
end


function frames = validateFrames(frames, nFrames)

    validateattributes(frames, {'numeric'}, ...
        {'vector', 'real', 'finite', 'integer', 'positive', 'nonempty'}, ...
        mfilename, 'framesToPlot', 2);

    frames = frames(:).';

    if any(frames > nFrames)
        error('saveKxyFramePNGs:FrameOutOfRange', ...
            'framesToPlot must contain values between 1 and %d.', nFrames);
    end
end


function [xLimits, yLimits] = getKxyLimits(kxy, frames)

    selected = kxy(:, :, frames);
    kx = real(selected);
    ky = imag(selected);

    finitePoints = isfinite(kx) & isfinite(ky);

    if ~any(finitePoints(:))
        error('saveKxyFramePNGs:NoFinitePoints', ...
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