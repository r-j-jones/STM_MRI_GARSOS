function outputs = plotGTimeAndFrequency(G, voxelA, voxelB, tPrime, dt)
%plotGTimeAndFrequency Plot g_{t'}(x,t) and g_{t'}(x,f).
%
% This function reproduces the g-only portions of Figure 4 in Lobos et al.
% It does not require the dynamic image series rho.
%
% Syntax:
%   plotGTimeAndFrequency(G, voxelA, voxelB)
%   plotGTimeAndFrequency(G, voxelA, voxelB, tPrime)
%   plotGTimeAndFrequency(G, voxelA, voxelB, tPrime, dt)
%   outputs = plotGTimeAndFrequency(...)
%
% Inputs:
%   G         Array of G(x) matrices with dimensions [Nx, Ny, T, T].
%             The final dimensions must be ordered as (tPrime,t), so that
%
%                 G(ix,iy,tPrime,t) = [G(x)]_{tPrime,t}.
%
%   voxelA    Spatial index [row, column] for voxel A.
%   voxelB    Spatial index [row, column] for voxel B.
%
% Optional inputs:
%   tPrime    Fixed t' used to extract g_{t'}(x,t). Default: min(75,T).
%
%   dt        Time between frames. Default: 1.
%             If dt=1, frequency is displayed in cycles/frame.
%             Otherwise, dt is assumed to be in seconds and frequency is
%             displayed in cycles/second.
%
% Output:
%   outputs   Struct containing:
%               time
%               frequency
%               gA
%               gB
%               gSpectrumA
%               gSpectrumB
%               voxelA
%               voxelB
%               tPrime
%               dt
%               figureHandle
%
% Mathematical relationship:
%
%   g_{t'}(x,t) = [G(x)]_{t',t}
%
%   g_{t'}(x,f) = FourierTransform_t{g_{t'}(x,t)}
%
% Example:
%   voxelA = [64, 38];
%   voxelB = [58, 52];
%   plotGTimeAndFrequency(G, voxelA, voxelB, 75, 1.6);
%
% Note:
%   If G is stored with the final dimensions ordered as (t,tPrime), modify
%   the two extraction statements as described near the relevant code.

    narginchk(3, 5);

    validateattributes(G, {'numeric'}, {'nonempty'});
    validateattributes(voxelA, {'numeric'}, ...
        {'vector', 'numel', 2, 'integer', 'positive', 'finite'});
    validateattributes(voxelB, {'numeric'}, ...
        {'vector', 'numel', 2, 'integer', 'positive', 'finite'});

    if ndims(G) ~= 4
        error('plotGTimeAndFrequency:InvalidGDimensions', ...
            'G must have dimensions [Nx, Ny, T, T].');
    end

    [Nx, Ny, T1, T2] = size(G);

    if T1 ~= T2
        error('plotGTimeAndFrequency:NonSquareG', ...
            ['The final two dimensions of G must be equal because each ' ...
             'G(x) matrix must have dimensions T-by-T.']);
    end

    T = T1;

    if nargin < 4 || isempty(tPrime)
        tPrime = min(75, T);
    end

    if nargin < 5 || isempty(dt)
        dt = 1;
    end

    validateattributes(tPrime, {'numeric'}, ...
        {'scalar', 'integer', '>=', 1, '<=', T, 'finite'});
    validateattributes(dt, {'numeric'}, ...
        {'scalar', 'real', 'positive', 'finite'});

    voxelA = voxelA(:).';
    voxelB = voxelB(:).';

    if voxelA(1) > Nx || voxelA(2) > Ny
        error('plotGTimeAndFrequency:VoxelAOutsideImage', ...
            'voxelA must be within [1,%d] by [1,%d].', Nx, Ny);
    end

    if voxelB(1) > Nx || voxelB(2) > Ny
        error('plotGTimeAndFrequency:VoxelBOutsideImage', ...
            'voxelB must be within [1,%d] by [1,%d].', Nx, Ny);
    end

    ixA = voxelA(1);
    iyA = voxelA(2);
    ixB = voxelB(1);
    iyB = voxelB(2);

    % Equation (17) in Lobos et al.:
    %
    %   g_{t'}(x,t) = [G(x)]_{t',t}.
    %
    % Hold t' fixed and extract row t' from each voxel's G(x) matrix.
    gA = reshape(G(ixA, iyA, tPrime, :), 1, T);
    gB = reshape(G(ixB, iyB, tPrime, :), 1, T);

    % If your G array is instead ordered as [Nx, Ny, t, tPrime], replace
    % the preceding two statements with:
    %
    % gA = reshape(G(ixA, iyA, :, tPrime), 1, T);
    % gB = reshape(G(ixB, iyB, :, tPrime), 1, T);

    % Temporal Fourier transforms. No spatial FFT is used here.
    gSpectrumA = fftshift(fft(gA, [], 2), 2);
    gSpectrumB = fftshift(fft(gB, [], 2), 2);

    % Time and temporal-frequency axes.
    time = (0:T-1) * dt;
    frequency = (-floor(T/2):ceil(T/2)-1) / (T * dt);

    if dt == 1
        timeLabel = 't [frames]';
        frequencyLabel = 'f [cycles/frame]';
    else
        timeLabel = 't [sec]';
        frequencyLabel = 'f [cycles/sec]';
    end

    % Use common limits so that voxel A and voxel B can be compared.
    timeMaximum = max([abs(gA), abs(gB)]);

    spectrumMaximum = max([ ...
        abs(gSpectrumA), abs(gSpectrumB)]);

    if isempty(timeMaximum) || timeMaximum == 0
        timeMaximum = 1;
    end

    if isempty(spectrumMaximum) || spectrumMaximum == 0
        spectrumMaximum = 1;
    end

    figureHandle = figure( ...
        'Color', 'w', ...
        'Name', 'Lobos et al. Figure 4: g-only plots');

    layout = tiledlayout(2, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    % Voxel A: time-domain magnitude.
    nexttile(layout, 1);
    plot(time, abs(gA), 'r-', 'LineWidth', 1.25);
    grid on;
    box on;
    xlim([time(1), time(end)]);
    ylim([0, 1.05 * timeMaximum]);
    xlabel(timeLabel);
    ylabel(sprintf('|g_{%d}(\\mathbf{x}_A,t)|', tPrime), ...
        'Interpreter', 'tex');
    title(sprintf('Voxel A: [%d, %d]', ixA, iyA));

    % Voxel B: time-domain magnitude.
    nexttile(layout, 2);
    plot(time, abs(gB), 'r-', 'LineWidth', 1.25);
    grid on;
    box on;
    xlim([time(1), time(end)]);
    ylim([0, 1.05 * timeMaximum]);
    xlabel(timeLabel);
    ylabel(sprintf('|g_{%d}(\\mathbf{x}_B,t)|', tPrime), ...
        'Interpreter', 'tex');
    title(sprintf('Voxel B: [%d, %d]', ixB, iyB));

    % Voxel A: temporal-frequency magnitude.
    nexttile(layout, 3);
    plot(frequency, abs(gSpectrumA), 'r-', 'LineWidth', 1.25);
    grid on;
    box on;
    xlim([frequency(1), frequency(end)]);
    ylim([0, 1.05 * spectrumMaximum]);
    xlabel(frequencyLabel);
    ylabel(sprintf('|\\breve{g}_{%d}(\\mathbf{x}_A,f)|', tPrime), ...
        'Interpreter', 'tex');
    title('Voxel A temporal spectrum');

    % Voxel B: temporal-frequency magnitude.
    nexttile(layout, 4);
    plot(frequency, abs(gSpectrumB), 'r-', 'LineWidth', 1.25);
    grid on;
    box on;
    xlim([frequency(1), frequency(end)]);
    ylim([0, 1.05 * spectrumMaximum]);
    xlabel(frequencyLabel);
    ylabel(sprintf('|\\breve{g}_{%d}(\\mathbf{x}_B,f)|', tPrime), ...
        'Interpreter', 'tex');
    title('Voxel B temporal spectrum');

    title(layout, sprintf( ...
        'Time and frequency behavior of g_{t''}(x,\\cdot), t'' = %d', ...
        tPrime), ...
        'Interpreter', 'tex');

    outputs = struct();
    outputs.time = time;
    outputs.frequency = frequency;
    outputs.gA = gA;
    outputs.gB = gB;
    outputs.gSpectrumA = gSpectrumA;
    outputs.gSpectrumB = gSpectrumB;
    outputs.voxelA = voxelA;
    outputs.voxelB = voxelB;
    outputs.tPrime = tPrime;
    outputs.dt = dt;
    outputs.figureHandle = figureHandle;
end
