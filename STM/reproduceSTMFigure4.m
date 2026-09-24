function reproduceSTMFigure4( rho, G, voxelA, voxelB, dt, tPrime )
%reproduceSTMFigure4 Reproduce the main plots in Lobos et al., Figure 4
%
% Required workspace variables:
%   rho       Dynamic image series, size [Nx, Ny, T].
%   G         G(x) matrices, size [Nx, Ny, T, T], where
%             G(ix,iy,tPrime,t) = [G(x)]_{tPrime,t}.
%   voxelA    [row, column] index of voxel A.
%   voxelB    [row, column] index of voxel B.
%
% Optional settings:
%   dt        Temporal sampling interval. Use dt = 1 for cycles/frame, or
%             the acquisition TR/frame interval for cycles/second.
%   tPrime    Fixed t' used to extract g_{t'}(x,t). Figure 4 uses t'=75.
%
% Example:
%   voxelA = [64, 38];
%   voxelB = [58, 52];
%   dt = 1;          % Frequency in cycles/frame
%   tPrime = 100;

if nargin<5 || isempty(dt)
    dt = 1;
end
if nargin<6 || isempty(tPrime)
    tPrime = 100;
end

assert(ndims(rho) == 3, ...
    'rho must have dimensions [Nx, Ny, T].');
assert(ndims(G) == 4, ...
    'G must have dimensions [Nx, Ny, T, T].');
assert(exist('voxelA', 'var') == 1 && numel(voxelA) == 2, ...
    'Define voxelA as [row, column].');
assert(exist('voxelB', 'var') == 1 && numel(voxelB) == 2, ...
    'Define voxelB as [row, column].');

[Nx, Ny, T] = size(rho);

assert(isequal(size(G), [Nx, Ny, T, T]), ...
    'G must have dimensions [Nx, Ny, T, T].');
assert(tPrime >= 1 && tPrime <= T && tPrime == fix(tPrime), ...
    'tPrime must be an integer between 1 and T.');
assert(all(voxelA >= 1) && voxelA(1) <= Nx && voxelA(2) <= Ny, ...
    'voxelA is outside the image.');
assert(all(voxelB >= 1) && voxelB(1) <= Nx && voxelB(2) <= Ny, ...
    'voxelB is outside the image.');

% Extract voxel indices.
ixA = voxelA(1);
iyA = voxelA(2);
ixB = voxelB(1);
iyB = voxelB(2);

% Extract the dynamic voxel signals.
rhoA = reshape(rho(ixA, iyA, :), 1, T);
rhoB = reshape(rho(ixB, iyB, :), 1, T);

% Equation (17):
%   g_{t'}(x,t) = [G(x)]_{t',t}.
%
% Hold t' fixed and extract ROW t' from each voxel's G(x) matrix.
gA = reshape(G(ixA, iyA, tPrime, :), 1, T);
gB = reshape(G(ixB, iyB, tPrime, :), 1, T);

% Calculate temporal spectra.
rhoSpectrumA = fftshift(fft(rhoA, [], 2), 2);
rhoSpectrumB = fftshift(fft(rhoB, [], 2), 2);
gSpectrumA   = fftshift(fft(gA,   [], 2), 2);
gSpectrumB   = fftshift(fft(gB,   [], 2), 2);

% Frequency axis corresponding to fftshift(fft(...)).
f = (-floor(T/2):ceil(T/2)-1) / (T * dt);
time = (0:T-1) * dt;

if dt == 1
    timeLabel = 't [frames]';
    frequencyLabel = 'f [cycles/frame]';
else
    timeLabel = 't [sec]';
    frequencyLabel = 'f [cycles/sec]';
end

% Display limits shared by corresponding A/B panels.
timeYMax = max([ ...
    abs(rhoA), abs(rhoB), abs(gA), abs(gB)], [], 'all');
spectrumYMax = max([ ...
    abs(rhoSpectrumA), abs(rhoSpectrumB), ...
    abs(gSpectrumA), abs(gSpectrumB)], [], 'all');

if timeYMax == 0
    timeYMax = 1;
end
if spectrumYMax == 0
    spectrumYMax = 1;
end

% Create a Figure 4-style layout.
figure('Color', 'w', 'Name', 'Lobos et al. STM Figure 4');
layout = tiledlayout(2, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

% Panel (a): first time frame and selected voxels.
nexttile(layout, 1, [2, 1]);
imagesc(abs(rho(:, :, 1)));
axis image;
axis off;
colormap(gray);
hold on;

plot(iyA, ixA, 'yo', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'y');
plot(iyB, ixB, 'yo', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'y');

text(iyA - 3, ixA, 'A', ...
    'Color', 'y', ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right');

text(iyB + 3, ixB, 'B', ...
    'Color', 'y', ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');

title('(a) Voxel locations');

% Panel (b), top: voxel A in time.
nexttile(layout, 2);
plot(time, abs(rhoA), 'b-', 'LineWidth', 1.1);
hold on;
plot(time, abs(gA), 'r-', 'LineWidth', 1.1);
grid on;
box on;
xlim([time(1), time(end)]);
ylim([0, 1.05 * timeYMax]);
xlabel(timeLabel);
title('Voxel A');
legend( ...
    {'|\rho(\mathbf{x}_A,t)|', ...
     sprintf('|g_{%d}(\\mathbf{x}_A,t)|', tPrime)}, ...
    'Location', 'best', ...
    'Interpreter', 'tex');

% Panel (c), top: voxel B in time.
nexttile(layout, 3);
plot(time, abs(rhoB), 'b-', 'LineWidth', 1.1);
hold on;
plot(time, abs(gB), 'r-', 'LineWidth', 1.1);
grid on;
box on;
xlim([time(1), time(end)]);
ylim([0, 1.05 * timeYMax]);
xlabel(timeLabel);
title('Voxel B');
legend( ...
    {'|\rho(\mathbf{x}_B,t)|', ...
     sprintf('|g_{%d}(\\mathbf{x}_B,t)|', tPrime)}, ...
    'Location', 'best', ...
    'Interpreter', 'tex');

% Panel (b), bottom: voxel A temporal spectra.
nexttile(layout, 5);
plot(f, abs(rhoSpectrumA), 'b-', 'LineWidth', 1.1);
hold on;
plot(f, abs(gSpectrumA), 'r-', 'LineWidth', 1.1);
grid on;
box on;
xlim([f(1), f(end)]);
ylim([0, 1.05 * spectrumYMax]);
xlabel(frequencyLabel);
legend( ...
    {'|\breve{\rho}(\mathbf{x}_A,f)|', ...
     sprintf('|\\breve{g}_{%d}(\\mathbf{x}_A,f)|', tPrime)}, ...
    'Location', 'best', ...
    'Interpreter', 'tex');

% Panel (c), bottom: voxel B temporal spectra.
nexttile(layout, 6);
plot(f, abs(rhoSpectrumB), 'b-', 'LineWidth', 1.1);
hold on;
plot(f, abs(gSpectrumB), 'r-', 'LineWidth', 1.1);
grid on;
box on;
xlim([f(1), f(end)]);
ylim([0, 1.05 * spectrumYMax]);
xlabel(frequencyLabel);
legend( ...
    {'|\breve{\rho}(\mathbf{x}_B,f)|', ...
     sprintf('|\\breve{g}_{%d}(\\mathbf{x}_B,f)|', tPrime)}, ...
    'Location', 'best', ...
    'Interpreter', 'tex');

title(layout, sprintf( ...
    'Temporal behavior and spectra for t'' = %d', tPrime));

% The extracted quantities are retained in the workspace:
%   gA, gB                         Time-domain g_{t'}(x,t)
%   rhoSpectrumA, rhoSpectrumB     Temporal spectra of rho(x,t)
%   gSpectrumA, gSpectrumB         Temporal spectra of g_{t'}(x,t)