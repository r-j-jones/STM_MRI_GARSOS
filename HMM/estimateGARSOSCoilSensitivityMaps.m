function [noiseWhiteningTransform, channelSensitivityMaps, adcChannelIds, debugStructure] = ...
    estimateGARSOSCoilSensitivityMaps(s)
%estimateGARSOSCoilSensitivityMaps Estimate receive-coil sensitivity/combination maps.
%
% This function processes a calibration measurement containing:
%   1. A surface-coil calibration acquisition
%   2. A body-coil calibration acquisition
%   3. Two surface-coil noise acquisitions
%   4. An end-of-scan packet
%
% The processing pipeline:
%   1. Validates the expected acquisition structure.
%   2. Separates surface-coil, body-coil, and noise readouts.
%   3. Assembles the surface- and body-coil calibration k-space arrays.
%   4. Estimates the surface-coil noise covariance matrix.
%   5. Constructs and applies a channel-noise whitening transform.
%   6. Reconstructs surface- and body-coil calibration images.
%   7. Estimates smooth, body-coil-referenced channel maps.
%   8. Assigns spatial metadata and generates diagnostic volumes.
%
%
% INPUT
%   s
%       Raw calibration data structure containing Siemens measurement data
%       headers, channel headers, raw samples, and measurement parameters.
%
% OUTPUTS
%   noiseWhiteningTransform
%       Matrix that whitens the surface-coil channel noise covariance.
%
%   channelSensitivityMaps
%       ArrayVolume containing the complex per-channel maps.
%
%   adcChannelIds
%       ADC channel identifiers for the surface-coil channels.
%
%   debugStructure
%       Structure containing body-coil and map-combination diagnostics.


%% Validate the expected calibration acquisition structure

% The calibration scan is expected to contain exactly two ICE acquisition
% indices. Acquisition 0 contains calibration and first-noise data, while
% acquisition 1 contains the second noise measurement.
nAcquisitions = numel(unique(s.m{1}.mdh.sLC.ushAcquisition));

assert(nAcquisitions == 2, ...
    'Expected two ICE acquisitions in calibration scan.');

% A zero-sample header at the end is used as the scan-termination packet.
assert(s.m{1}.mdh.ushSamplesInScan(end) == 0, ...
    'Calibration scan terminated unexpectedly (incomplete data).');

% Acquisition 0 must contain two ICE sets:
%   set 0: surface-coil calibration
%   set 1: body-coil calibration
nSets = numel(unique( ...
    s.m{1}.mdh.sLC.ushSet( ...
    s.m{1}.mdh.sLC.ushAcquisition == 0)));

assert(nSets == 2, ...
    'Expected two ICE sets in calibration scan.');


%% Identify surface-coil, body-coil, noise, and end-of-scan packets

% Surface-coil calibration readouts have 128 samples and belong to
% acquisition 0, set 0.
ii1 = find( ...
    s.m{1}.mdh.sLC.ushAcquisition == 0 & ...
    s.m{1}.mdh.sLC.ushSet == 0 & ...
    s.m{1}.mdh.ushSamplesInScan == 128);

% Body-coil calibration readouts have 128 samples and belong to
% acquisition 0, set 1.
ii2 = find( ...
    s.m{1}.mdh.sLC.ushAcquisition == 0 & ...
    s.m{1}.mdh.sLC.ushSet == 1 & ...
    s.m{1}.mdh.ushSamplesInScan == 128);

% First noise measurement: acquisition 0, set 0, with 512 samples.
ii3 = find( ...
    s.m{1}.mdh.sLC.ushAcquisition == 0 & ...
    s.m{1}.mdh.sLC.ushSet == 0 & ...
    s.m{1}.mdh.ushSamplesInScan == 512);

% Second noise measurement: acquisition 1, set 0, with 512 samples.
ii4 = find( ...
    s.m{1}.mdh.sLC.ushAcquisition == 1 & ...
    s.m{1}.mdh.sLC.ushSet == 0 & ...
    s.m{1}.mdh.ushSamplesInScan == 512);

% End-of-scan packet.
ii5 = find(s.m{1}.mdh.ushSamplesInScan == 0);

% Count each readout category.
n1 = numel(ii1);
n2 = numel(ii2);
n3 = numel(ii3);
n4 = numel(ii4);
n5 = numel(ii5);

% Surface- and body-coil calibration scans must contain equal numbers of
% readouts so that they represent matching calibration grids.
assert(n1 == n2);

% The two noise measurements must contain equal numbers of readouts.
assert(n3 == n4);

% Confirm that every raw-data header belongs to exactly one of the expected
% categories above. This rejects unexpected packet or readout types.
nReadouts = size(s.m{1}.mdh, 1);

assert(n1 + n2 + n3 + n4 + n5 == nReadouts);


%% Extract calibration-grid coordinates

% The scanner's line and partition counters locate each surface-coil
% readout in the 32-by-32 calibration grid.
x1 = s.m{1}.mdh.sLC.ushLine(ii1);
y1 = s.m{1}.mdh.sLC.ushPartition(ii1);

% Corresponding locations for the body-coil calibration readouts.
x2 = s.m{1}.mdh.sLC.ushLine(ii2);
y2 = s.m{1}.mdh.sLC.ushPartition(ii2);


%% Determine the number of receive channels

% Each channelHeaders entry contains one row per acquired receive channel.
nSurfaceCoils = size( ...
    s.m{1}.mdh.channelHeaders{ii1(1)}, 1);

nBodyCoils = size( ...
    s.m{1}.mdh.channelHeaders{ii2(1)}, 1);


%% Assemble surface-coil calibration k-space

% Convert zero-based scanner line/partition counters to MATLAB's one-based
% linear indices for a 32-by-32 calibration grid.
jj1 = sub2ind([32 32], x1 + 1, y1 + 1);

% Each populated cell will contain one 128-sample readout across all
% surface-coil channels.
c1 = cell(32, 32);

% Convert each readout's raw channel data into the dimension ordering
% expected by the subsequent cell-array concatenation.
c1(jj1) = cellfun( ...
    @(x) permute(cell2mat(x.rawdata'), [3 1 4 2]), ...
    s.m{1}.mdh.channelHeaders(ii1), ...
    'UniformOutput', false);

% Record which calibration-grid locations contain measured data.
m1 = ~cellfun(@(x) isempty(x), c1);

% Explicitly zero-fill missing calibration-grid locations. Each cell has
% dimensions compatible with one readout and all surface-coil channels.
c1(~m1) = {zeros(1, 128, 1, nSurfaceCoils, 'single')};

% Rearrange and concatenate the cells into a numerical k-space array. The
% fourth dimension is the surface-coil channel dimension.
C1 = cell2mat(permute(c1, [1 3 2 4]));


%% Assemble body-coil calibration k-space

% Convert body-coil line/partition counters to linear grid indices.
jj2 = sub2ind([32 32], x2 + 1, y2 + 1);

% Allocate the body-coil calibration grid.
c2 = cell(32, 32);

% Insert each acquired body-coil readout into its calibration-grid location.
c2(jj2) = cellfun( ...
    @(x) permute(cell2mat(x.rawdata'), [3 1 4 2]), ...
    s.m{1}.mdh.channelHeaders(ii2), ...
    'UniformOutput', false);

% Identify measured grid locations.
m2 = ~cellfun(@(x) isempty(x), c2);

% Zero-fill unmeasured locations using the body-coil channel count.
c2(~m2) = {zeros(1, 128, 1, nBodyCoils, 'single')};

% Concatenate the cell array into a numerical body-coil k-space array.
C2 = cell2mat(permute(c2, [1 3 2 4]));


%% Assemble the two surface-coil noise measurements

% Convert the first noise acquisition into a numerical array and reshape it
% so that each row is one noise sample and each column is one surface-coil
% channel.
c3 = cell2mat(cellfun( ...
    @(x) permute(cell2mat(x.rawdata'), [3 1 4 2]), ...
    s.m{1}.mdh.channelHeaders(ii3), ...
    'UniformOutput', false));

C3 = reshape(c3, [128 * 512, nSurfaceCoils]);

% Process the second noise acquisition using the same channel ordering.
c4 = cell2mat(cellfun( ...
    @(x) permute(cell2mat(x.rawdata'), [3 1 4 2]), ...
    s.m{1}.mdh.channelHeaders(ii4), ...
    'UniformOutput', false));

C4 = reshape(c4, [128 * 512, nSurfaceCoils]);

% Combine both noise measurements to increase the number of samples used
% for estimating the inter-channel covariance.
C34 = cat(1, C3, C4);


%% Estimate the noise covariance and construct the whitening transform

% Estimate the complex surface-coil noise covariance matrix. Conjugate
% transpose is supplied by MATLAB's apostrophe operator.
noiseCovarianceMatrix = C34' * C34 / size(C34, 1);

% MATLAB's chol returns an upper-triangular matrix R satisfying:
%   noiseCovarianceMatrix = R' * R
noiseColoringTransform = chol(noiseCovarianceMatrix);

% The inverse coloring transform whitens the channel noise when applied
% across the surface-coil channel dimension.
noiseWhiteningTransform = inv(noiseColoringTransform);


%% Apply noise whitening to the surface-coil calibration data

% Contract the fourth dimension of C1 with the first dimension of the
% whitening matrix. The output's fourth dimension contains the whitened
% channel basis.
C1decorr = tprod( ...
    C1, [1 2 3 -5], ...
    noiseWhiteningTransform, [-5 4]);

% Release the original unwhitened surface-coil calibration array.
clear C1;


%% Construct the calibration k-space Gaussian window

% Gaussian width expressed in cycles per field of view.
sigma1 = 15 / 2;

% Coordinates match the 32-by-128-by-32 calibration array. The second
% coordinate is scaled by one-half to account for its different sampling
% extent/convention.
[I, J, K] = ndgrid( ...
    -16:15, ...
    (-64:63) / 2, ...
    -16:15);

% Three-dimensional Gaussian weighting function for calibration k-space.
F1 = exp(-(I.^2 + J.^2 + K.^2) / (2 * sigma1^2));


%% Define reconstruction array sizes

% Obtain the nominal scanner base resolution.
baseResolution = double( ...
    s.m{1}.param.Meas{1}.MEAS.sKSpace.lBaseResolution);

% Nominal map dimensions before readout oversampling and zero-padding.
nominalMapArraySize = [ ...
    baseResolution, ...
    baseResolution, ...
    baseResolution];

% Preserve a factor-of-two size in the second dimension, corresponding to
% the reconstruction's readout oversampling convention.
reconstructionMapArraySize = [ ...
    baseResolution, ...
    2 * baseResolution, ...
    baseResolution];

% Double all three spatial dimensions by centered resizing/zero-padding.
% Preserve the appropriate channel count as dimension four.
newArraySize1 = [ ...
    2 * reconstructionMapArraySize, ...
    size(C1decorr, 4)];

newArraySize2 = [ ...
    2 * reconstructionMapArraySize, ...
    size(C2, 4)];


%% Construct the Gaussian used to smooth the map estimator

sigma2 = 15 / 2;

% Coordinates span the final surface-coil reconstruction dimensions. The
% second coordinate again uses a factor-of-two coordinate scaling.
[I, J, K] = ndgrid( ...
    -(newArraySize1(1) / 2):(newArraySize1(1) / 2 - 1), ...
    (-(newArraySize1(2) / 2):(newArraySize1(2) / 2 - 1)) / 2, ...
    -(newArraySize1(3) / 2):(newArraySize1(3) / 2 - 1));

% Frequency-domain Gaussian used below to smooth the numerator and
% denominator of the channel-map estimator.
F2 = exp(-(I.^2 + J.^2 + K.^2) / (2 * sigma2^2));


%% Reconstruct reference calibration images

% Centered resizing zero-fills or crops the calibration k-space to the
% requested output size. The centered FFT transforms the first three
% dimensions according to this codebase's reconstruction convention.
D1ref = centeredFFT( ...
    centeredResize( ...
    bsxfun(@times, C1decorr, 1), ...
    newArraySize1), ...
    3);

D2ref = centeredFFT( ...
    centeredResize( ...
    bsxfun(@times, C2, 1), ...
    newArraySize2), ...
    3);


%% Reconstruct the calibration images used for map estimation

% % Reconstruct a Gaussian-windowed surface-coil calibration image.
% D1 = centeredFFT( ...
%     centeredResize( ...
%     bsxfun(@times, C1decorr, F1), ...
%     newArraySize1), ...
%     3);

% Preserve the original implementation exactly: the preceding D1 result is
% immediately replaced by a reconstruction without F1 weighting. Thus, the
% surface-coil image ultimately used below is not weighted by F1.
D1 = centeredFFT( ...
    centeredResize( ...
    bsxfun(@times, C1decorr, 1), ...
    newArraySize1), ...
    3);

% The body-coil image used for map estimation retains the F1 Gaussian
% calibration-k-space weighting.
D2 = centeredFFT( ...
    centeredResize( ...
    bsxfun(@times, C2, F1), ...
    newArraySize2), ...
    3);


%% Combine the body-coil channels

% Coherently sum all body-coil channels. The estimation below uses the
% magnitude of D20 as the relatively homogeneous reference image.
D20 = sum(D2, 4);

% Also form an unfiltered/interpolated body-coil reference for diagnostics.
D20ref = sum(D2ref, 4);


%% Form and smooth the channel-map estimator

% For each whitened surface-coil channel, form a body-reference-weighted
% conjugate image. This is the numerator of a normalized local
% coil-combination estimate.
J0 = bsxfun(@times, conj(D1), abs(D20));

% Sum the whitened surface-coil image power across channels. This common
% denominator normalizes the per-channel conjugate weights.
W0 = sum(abs(D1).^2, 4);

% Smooth the estimator numerator by transforming, multiplying by the
% Gaussian F2, and transforming back. The first three dimensions are the
% spatial dimensions; the fourth channel dimension is preserved.
JF0 = centeredFFT( ...
    bsxfun(@times, centeredIFFT(J0, 3), F2), ...
    3);

% Apply the same smoothing operation to the common power denominator.
WF0 = centeredFFT( ...
    bsxfun(@times, centeredIFFT(W0, 3), F2), ...
    3);

% Divide each channel's smoothed numerator by the smoothed total channel
% power. No additional mask, threshold, or denominator regularization is
% introduced, preserving the original operation.
C = bsxfun(@rdivide, JF0, WF0);


%% Generate a diagnostic surface-coil combination

% Apply the estimated maps to the unfiltered surface-coil reference images
% and sum across channels. This image should resemble the body-coil
% reference magnitude when the maps behave as intended.
T2 = sum(C .* D1ref, 4);


%% Record the original surface-coil ADC channel identifiers

adcChannelIds = double( ...
    s.m{1}.mdh.channelHeaders{ii1(1)}.ushChannelId);


%% Store the maps in an ArrayVolume

channelSensitivityMaps = ArrayVolume(C);


%% Assign spatial field-of-view and voxel-size metadata

% Scanner field-of-view values are ordered here as phase, readout, and
% slice/partition.
nominalFieldOfViewPRS = [ ...
    s.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dPhaseFOV, ...
    s.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dReadoutFOV, ...
    s.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dThickness];

% The factor of two reflects the doubled map dimensions produced by the
% centered resizing above.
voxelSize = ( ...
    nominalFieldOfViewPRS ./ nominalMapArraySize / 2 ...
    )';

% Preserve the original hard-coded mapping from array coordinates to the
% spatial coordinate system used by ArrayVolume.
channelSensitivityMaps.R = [ ...
     0  0  1; ...
    -1  0  0; ...
     0 -1  0];

% Assign voxel spacing.
channelSensitivityMaps.v = voxelSize;

% Set the spatial origin using the reconstruction-array center and the
% original one-voxel adjustment in the third coordinate.
channelSensitivityMaps.r0 = ...
    channelSensitivityMaps.R * ...
    (diag(voxelSize) * ...
    (-floor(2 * reconstructionMapArraySize' / 2) + [0 0 -1]'));


%% Create diagnostic ArrayVolume objects with matching geometry

% Complex, coherently combined body-coil reference.
debugStructure.D20ref = ...
    channelSensitivityMaps.cloneWithNewContent(D20ref);

% Surface-coil image formed using the estimated channel maps.
debugStructure.T2 = ...
    channelSensitivityMaps.cloneWithNewContent(T2);

% Magnitude of the body-coil reference.
debugStructure.bodyCoilImage = ...
    channelSensitivityMaps.cloneWithNewContent(abs(D20ref));

% debugStructure.C1 = C1;
% debugStructure.C1decorr = C1decorr;

debugStructure.nSurfaceCoils = nSurfaceCoils;
debugStructure.nBodyCoils = nBodyCoils;

debugStructure.noiseCovarianceMatrix = noiseCovarianceMatrix;
debugStructure.noiseWhiteningTransform = noiseWhiteningTransform;
end