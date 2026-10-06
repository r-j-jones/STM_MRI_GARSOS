function figuresOut = plotGARSOSKspaceCoverage(results, varargin)
%plotGARSOSKSPACE COVERAGE Plot trajectory-derived coverage diagnostics.

    p = inputParser;
    addParameter(p, 'OverlaySamples', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    figuresOut = gobjects(0);
    figuresOut(end + 1) = figure('Name', 'GARSOS raw k-space coverage');
    imagesc(results.raw.normalizedMap);
    axis image;
    colorbar;
    title(sprintf('Raw coverage H / mean(H), CV = %.2f%%', ...
        results.raw.CV_percent));
    xlabel('kx bin');
    ylabel('ky bin');
    if p.Results.OverlaySamples
        hold on;
        plot(results.sampleCoordinates(:, 1), results.sampleCoordinates(:, 2), ...
            '.', 'Color', [1 1 1], 'MarkerSize', 3);
        hold off;
    end

    figuresOut(end + 1) = figure('Name', 'GARSOS DCF-weighted coverage');
    imagesc(results.dcf.normalizedMap);
    axis image;
    colorbar;
    title(sprintf('DCF-weighted coverage Hw / mean(Hw), CV = %.2f%%', ...
        results.dcf.CV_percent));
    xlabel('kx bin');
    ylabel('ky bin');

    figuresOut(end + 1) = figure('Name', 'GARSOS nearest-neighbor distances');
    histogram(results.nearestNeighbor.distances);
    xlabel('Nearest-neighbor distance');
    ylabel('Samples');
    title('Trajectory nearest-neighbor distances');

    figuresOut(end + 1) = figure('Name', 'GARSOS angular gaps');
    hold on;
    for iBand = 1:numel(results.angular)
        if ~isempty(results.angular(iBand).gaps)
            histogram(results.angular(iBand).gaps, ...
                'DisplayStyle', 'stairs', ...
                'DisplayName', sprintf('%.0f-%.0f%% radius', ...
                100 * results.angular(iBand).radialRange));
        end
    end
    hold off;
    xlabel('Angular gap (radians)');
    ylabel('Count');
    title('Angular-gap distributions');
    legend('show');
end
