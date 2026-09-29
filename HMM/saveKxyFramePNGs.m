function saveKxyFramePNGs(kxy, framesToPlot, outputDirectory, zoomLim)
%SAVEKXYFRAMEPNGS Save full-view and zoomed radial trajectory PNGs.
%
%   saveKxyFramePNGs(kxy, framesToPlot, outputDirectory, zoomLim)
%
% Inputs
%   kxy
%       Complex numeric array of size:
%       nReadout x nSpokesPerFrame x nFrames
%
%       real(kxy) contains kx and imag(kxy) contains ky.
%
%   framesToPlot
%       Vector of frame indices to plot.
%
%   outputDirectory
%       Directory in which the PNG files are saved.
%
%   zoomLim
%       Positive scalar defining the zoomed limits:
%       xlim([-zoomLim zoomLim])
%       ylim([-zoomLim zoomLim])
%
% Example
%   saveKxyFramePNGs(kxy, 1:133, 'trajectory_pngs', 0.05);
%
% USAGE:
%  framesToPlot = 1:133;
%  outputDirectory = fullfile(pwd, 'trajectory_plots');
%  zoomLim = 0.05;
%  saveKxyFramePNGs(kxy, framesToPlot, outputDirectory, zoomLim);

    narginchk(4, 4);

    validateattributes(kxy, {'numeric'}, ...
        {'nonempty'}, mfilename, 'kxy', 1);

    if ndims(kxy) ~= 3
        error('saveKxyFramePNGs:InvalidKxyDimensions', ...
            ['kxy must have dimensions nReadout x ' ...
             'nSpokesPerFrame x nFrames.']);
    end

    nSpokesPerFrame = size(kxy, 2);
    nFrames = size(kxy, 3);

    validateattributes(framesToPlot, {'numeric'}, ...
        {'vector', 'real', 'finite', 'integer', 'positive', 'nonempty'}, ...
        mfilename, 'framesToPlot', 2);

    framesToPlot = framesToPlot(:).';

    if any(framesToPlot > nFrames)
        error('saveKxyFramePNGs:FrameOutOfRange', ...
            'Frame indices must be between 1 and %d.', nFrames);
    end

    validateattributes(zoomLim, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'zoomLim', 4);

    if ~(ischar(outputDirectory) || ...
            (isstring(outputDirectory) && isscalar(outputDirectory)))
        error('saveKxyFramePNGs:InvalidOutputDirectory', ...
            'outputDirectory must be a character vector or string scalar.');
    end

    outputDirectory = char(outputDirectory);

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    % Use common full-view limits for all requested frames.
    selectedKxy = kxy(:, :, framesToPlot);
    selectedKx = real(selectedKxy);
    selectedKy = imag(selectedKxy);

    finitePoints = isfinite(selectedKx) & isfinite(selectedKy);

    if ~any(finitePoints(:))
        error('saveKxyFramePNGs:NoFinitePoints', ...
            'The requested frames contain no finite trajectory points.');
    end

    selectedKx = selectedKx(finitePoints);
    selectedKy = selectedKy(finitePoints);

    commonMinimum = min([selectedKx(:); selectedKy(:)]);
    commonMaximum = max([selectedKx(:); selectedKy(:)]);
    commonRange = commonMaximum - commonMinimum;

    if commonRange == 0
        padding = max(1, abs(commonMinimum)) * 0.05;
    else
        padding = 0.05 * commonRange;
    end

    % Identical x and y limits ensure a square plotting region.
    fullLimits = [commonMinimum - padding, commonMaximum + padding];

    fig = figure( ...
        'Color', 'w', ...
        'Visible', 'off', ...
        'Position', [2045 296 496 413]); %[1265 295 496 413]);

    colors = lines(nSpokesPerFrame);

    for iFrame = framesToPlot
        clf(fig);

        ax = axes(fig);
        hold(ax, 'on');

        % Initially draw lines without point markers.
        for iSpoke = 1:nSpokesPerFrame
            trajectory = kxy(:, iSpoke, iFrame);

            plot(ax, real(trajectory), imag(trajectory), ...
                '-', ...
                'Color', 'k', ...
                'LineWidth', 1.5);
        end

        xlabel(ax, 'k_x');
        ylabel(ax, 'k_y');
        title(ax, sprintf('Frame %d', iFrame));
        set(ax,'FontSize',15);

        xlim(ax, fullLimits);
        ylim(ax, fullLimits);
        axis(ax, 'equal');
        axis(ax, 'square');
        grid(ax, 'on');
        box(ax, 'on');

        % Save the original full-view PNG.
        fullOutputPath = fullfile(outputDirectory, ...
            sprintf('kxy_frame_%04d.png', iFrame));

        figure(fig); drawnow;
        print(gcf, fullOutputPath, '-dpng', '-r300');

        % Add a circular marker at every trajectory sample.
        cla(ax);
        for iSpoke = 1:nSpokesPerFrame
            trajectory = kxy(:, iSpoke, iFrame);

            plot(ax, real(trajectory), imag(trajectory), ...
                'o', ...
                'LineStyle', 'none', ...
                'Color', 'k', ...
                'MarkerSize', 4, ...
                'LineWidth', 0.5);
        end

        % Apply equal, square central zoom limits.
        axis(ax, 'equal');
        axis(ax, 'square');
        xlim(ax, [-zoomLim, zoomLim]);
        ylim(ax, [-zoomLim, zoomLim]);

        title(ax, sprintf('Frame %d zoom-in', iFrame));
        set(ax,'FontSize',15);

        % Save the zoomed version separately.
        zoomOutputPath = fullfile(outputDirectory, ...
            sprintf('kxy_frame_%04d_zoom.png', iFrame));

        figure(fig); drawnow;
        print(gcf, zoomOutputPath, '-dpng', '-r300');
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