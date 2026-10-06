I’ll trace the prepared 58-partition arrays through trajectory construction, gridding, FFT, resizing, and coil accumulation, then compare that exact path with the custom slice-wise IFFT.

Read , lines 1 to 150
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/GpuReconstructur20181008.m",
  "view_range": [
    1,
    150
  ]
}

Read , lines 1 to 180
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/shiftAndResize.m",
  "view_range": [
    1,
    180
  ]
}

Read , lines 175 to 275
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/reconstructGARFrameForSTM.m",
  "view_range": [
    175,
    275
  ]
}

Read , lines 55 to 140
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/loadPrepareGARSOSData.m",
  "view_range": [
    55,
    140
  ]
}

Read , lines 35 to 145
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/prepareRadialVibeDataForGridding2_3.m",
  "view_range": [
    35,
    145
  ]
}

Read , lines 65 to 165
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/prepareGARSOSSliceDataForMCNUFFT.m",
  "view_range": [
    65,
    165
  ]
}

Read , lines 1 to 160
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/centeredArrayVolumeResize.m",
  "view_range": [
    1,
    160
  ]
}

Read , lines 1 to 180
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/resampleArrayVolumeToSame.m",
  "view_range": [
    1,
    180
  ]
}

Read , lines 1 to 180
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/apodization3d2.m",
  "view_range": [
    1,
    180
  ]
}

Search for `imageSize|gridOversamplingFactor|nPartitions|nFourierPartitions|gridSize`

Read , lines 1 to 140
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/centeredResize.m",
  "view_range": [
    1,
    140
  ]
}

Read , lines 1 to 120
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/prepareGARSOSCalibrationData.m",
  "view_range": [
    1,
    120
  ]
}

Read , lines 1 to 80
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/prepareGARSOSSliceDataForMCNUFFT.m",
  "view_range": [
    1,
    80
  ]
}

Read , lines 1 to 180
Completed with input: {
  "path": "/Users/robertjones/Desktop/hmm_code_matlab/STM_MRI_GARSOS/HMM/validateGARSOSTrajectoryData.m",
  "view_range": [
    1,
    180
  ]
}

Search for `58|96|GpuReconstructur|coil|partition|slice|shiftAndResize`

# 1. Executive conclusion

## Exact 58 → 96 mechanism

The working reconstruction does **not** perform:

```text
58 → zero-fill to 96 → IFFT
```

It performs:

```text
58 nonuniform kz trajectory locations
    → 3-D Kaiser–Bessel gridding onto gridSize = [308 308 132]
    → 3-D unnormalized ifftn on the 308 × 308 × 132 grid
    → centered extraction/resizing to imageSize = [224 224 96]
```

The exact output-size transition occurs in `GpuReconstructur20181008.m`, where the gridded array has size `gridSize`, and then in:

```matlab
A = shiftAndResize(A, [], obj.imageSize, ...
    repmat('e', [1, obj.nDimensions]), ...
    repmat('m', [1, obj.nDimensions]));
```

with `obj.imageSize = [224 224 96]`.

For the stated parameters:

```matlab
imageSize = [224 224 96];
gridSize = ceil(1.375 .* imageSize);
```

therefore:

```text
gridSize = [308 308 132]
```

The z dimension is thus:

```text
58 supplied trajectory locations
→ 132-point Cartesian gridding grid
→ 132-point inverse FFT
→ 96-point centered output
```

The 58 locations are not treated as a regular 58-point image-space z array.

## Can the existing GPU backbone directly provide STM input?

Not directly.

`GpuReconstructur20181008`:

- accepts samples with all three trajectory coordinates;
- grids all three spatial dimensions;
- applies coil-combination profiles;
- sums all channels into one volume.

Its current output is therefore:

```text
224 × 224 × 96
```

not:

```text
224 × 224 × 96 × nCoils
```

It can be parameterized or minimally modified to retain per-coil gridded volumes, but those would be Cartesian image volumes. They would not be the radial STM input expected by `STM_reconstruction_GARSOS`, which requires data organized approximately as:

```text
nReadout × nSpokesPerFrame × nCoils
```

for a selected physical z slice.

The smallest reliable design is to reuse the **z-direction gridding/reconstruction logic**, while preserving the readout, spoke, and coil dimensions:

```text
[nReadout × 58 × nSpokes × nCoils]
    → z-direction nonuniform gridding to the 132-point Cartesian grid
    → z-direction IFFT
    → centered extraction to 96 z positions
    → select one z plane
    → [nReadout × nSpokesPerFrame × nCoils]
```

That preserves the in-plane radial coordinates and channel dimension. Reusing the full 3-D reconstructor and then taking a 2-D FFT of the resulting slice does not preserve the original radial sample locations.

---

# 2. Dimension and operation trace

## 2.1 Prepared data dimensions

The prepared raw data are assembled in `prepareRadialVibeDataForGridding2_3.m`:

```matlab
Z = decorrelateCoils2(Z, noiseWhiteningTransform);
Z = reshape(Z, [nLines*nPartitions*nSamples, nTransformedChannels]);
sampledData = double(Z);
```

The flattened `sampledData` layout is:

```text
[nReadout × nPartition × nSpokes, nCoils]
```

The user-reported prepared partition count is 58.

The trajectory is stored as:

```matlab
kSpaceLocations = [kx ky kz];
```

with each row corresponding to one flattened sample. The `kz` expression is constructed earlier in the same file:

```matlab
kz = (double(iiPartition(:)) - cPartition) ...
    / s.m{end}.param.Meas{1}.MEAS.sKSpace.lPartitions ...
    / s.m{end}.param.Meas{1}.MEAS.sKSpace.dSliceResolution;
```

at `prepareRadialVibeDataForGridding2_3.m`.

The code then repeats `kz` across readout samples and flattens it alongside `kx` and `ky`.

The earlier preparation also contains:

```matlab
Z = padarray(Z,[0 nFourierPartitions 0 0],0,'pre');
```

at `prepareRadialVibeDataForGridding2_3.m`.

This is zero-padding along the partition-related raw-data axis before flattening. The corresponding trajectory array is generated with the prepared partition layout and is passed unchanged into the later reconstruction. The important point for this review is that the reconstructor receives the prepared 58-partition data and their existing `kz` coordinates; it does not infer a regular 58-point image-space z grid.

## 2.2 Per-frame input to the 3-D reconstructor

In `reconstructGARFrameForSTM.m`, one temporal frame is selected:

```matlab
frameData = sampledData(sampleMask, :);
frameCoordinates = sampleCoordinates(sampleMask, :);
```

The resulting dimensions are approximately:

```text
frameData        [nReadout × nPartition × nSpokesPerFrame, nCoils]
frameCoordinates [nReadout × nPartition × nSpokesPerFrame, 3]
```

The coordinates are transformed:

```matlab
frameTrajectory = frameCoordinates * coordinateTransform;
```

at line 201.

Thus the reconstructor receives the transformed `frameTrajectory`, not necessarily the original untransformed `kx`, `ky`, and `kz`.

The code applies phase and density weighting:

```matlab
phaseCorrection = exp(1i * 2 * pi * ...
    (frameCoordinates * translationVector));

densityWeight = frameDCF * prod(visibleObjectSize(1:2)) ...
    * 2 / p.Results.PhaseWidth;

weightedFrameData = bsxfun(@times, frameData, ...
    densityWeight .* gaussianWeight .* phaseCorrection);
```

at `reconstructGARFrameForSTM.m`.

These operations change sample values, but not the number or organization of the trajectory locations.

---

# 3. Inside `GpuReconstructur20181008`

## 3.1 Grid dimensions

The constructor receives:

```matlab
gridSize
imageSize
kernelSize
sampleCoordinates
coilCombinationProfiles
```

at `GpuReconstructur20181008.m`.

For the reported configuration:

```matlab
imageSize = [224 224 96];
gridSize = [308 308 132];
kernelSize = [7 7 7];
```

The `gridSize` is constructed in `loadPrepareGARSOSData.m`:

```matlab
gridSize = ceil(griddingParams.gridOversamplingFactor .* imageSize);
```

with default oversampling factor `1.375`.

Therefore:

```text
input/output image grid: 224 × 224 × 96
oversampled Cartesian grid: 308 × 308 × 132
```

## 3.2 Mapping trajectory locations to the Cartesian grid

The constructor computes:

```matlab
kGrid = sampleCoordinates .* gridSize;
ii = round(kGrid);
dk = ii - kGrid;
```

at `GpuReconstructur20181008.m`.

For each sample, `kGrid` contains the three-dimensional Cartesian-grid coordinate corresponding to the transformed trajectory point.

The integer grid locations are then expanded by the kernel neighborhood:

```matlab
II = mod(reshape(permute(ii, [1 3 2]) + ...
    permute(pp, [3 2 1]), ...
    [nSamples*prod(kernelSize), nDimensions]), ...
    gridSize) * gridStride.' + 1;
```

at lines 56–56.

The `mod(..., gridSize)` operation wraps indices periodically. The resulting `II` contains the Cartesian grid bins receiving contributions from each sample and its `[7 7 7]` neighborhood.

The interpolation weights are computed independently in each dimension using the Kaiser–Bessel kernel:

```matlab
DK = ...
    interp1(kKernel{iDimension}, ...
    kernelSamples2{iDimension}, ...
    dk(:, iDimension) + p{iDimension});
```

and multiplied across dimensions.

Therefore, the 58 kz values are not assigned one-to-one to 58 output z planes. Each sample contributes to several z-grid locations according to its fractional grid position and the z component of the 3-D kernel.

## 3.3 Zero filling

Inside `grid`, the gridded Cartesian array is initialized indirectly by:

```matlab
A = accumarray(obj.II, DKi, [prod(obj.gridSize), 1]);
A = reshape(A, obj.gridSize);
```

at `GpuReconstructur20181008.m`.

`accumarray` fills Cartesian grid bins that receive no samples with zero. Thus there is zero-filled support on the 308 × 308 × 132 grid, but this is not equivalent to appending zeros to a 58-point regular kz vector.

The earlier `padarray` in `prepareRadialVibeDataForGridding2_3.m` also inserts zero-valued prepared partition entries, but the working reconstructor still uses the resulting trajectory coordinates and performs nonuniform 3-D gridding.

## 3.4 Inverse Fourier transform

After Cartesian deposition:

```matlab
A = ifftn(A);
```

at `GpuReconstructur20181008.m`.

This is a 3-D inverse FFT over:

```text
308 × 308 × 132
```

It is not an inverse FFT over the 58 partition entries.

The exact FFT normalization is MATLAB’s standard `ifftn` convention, with the usual inverse-transform scale over all three dimensions.

## 3.5 Conversion from 308 × 308 × 132 to 224 × 224 × 96

After `ifftn`, the code calls:

```matlab
A = shiftAndResize(A, [], obj.imageSize, ...
    repmat('e', [1, obj.nDimensions]), ...
    repmat('m', [1, obj.nDimensions]));
```

at line 105.

`shiftAndResize` computes:

```matlab
copySize = min([copySize; inputSize; outputSize], [], 1);
```

in `shiftAndResize.m`.

With:

```text
inputSize  = [308 308 132]
outputSize = [224 224 96]
```

the copy size is:

```text
copySize = [224 224 96]
```

Because the requested output is smaller than the input, the function extracts a centered/wrapped portion of the oversampled image grid. It uses:

```matlab
inputZeroLocations  = 'e'
outputZeroLocations = 'm'
```

where:

- `'e'` means the input zero index is at the first array element;
- `'m'` means the output zero index is centered.

The output is therefore a centered 224 × 224 × 96 image grid. There is no additional FFT at this stage and no interpolation from 58 to 96 image planes.

The exact size transition is therefore:

```text
308 × 308 × 132 after gridding and IFFT
→ 224 × 224 × 96 after shiftAndResize
```

## 3.6 Coil-combination boundary

The constructor computes:

```matlab
apodizationFunction = coilCombinationProfiles * prod(gridSize);
```

at `GpuReconstructur20181008.m`.

For each coil:

```matlab
DKi = reshape(obj.DK .* sampledData(:, iChannel), ...);
A = accumarray(...);
A = reshape(A, obj.gridSize);
A = ifftn(A);
A = shiftAndResize(...);
A = A(:);
A = A .* obj.apodizationFunction(:, iChannel);
```

at lines 93–112.

Then the channels are summed:

```matlab
if iChannel == 1
    imageData = A;
else
    imageData = imageData + A;
end
```

at lines 113–119.

Finally:

```matlab
imageData = reshape(imageData, obj.imageSize);
```

at line 123.

Thus the current output is:

```text
224 × 224 × 96
```

not:

```text
224 × 224 × 96 × nCoils
```

The coil-combination boundary is the accumulation:

```matlab
imageData = imageData + A;
```

The map/profile weighting happens before that accumulation.

---

# 4. Physical z mapping

## 4.1 Target geometry

The target image size is selected in `loadPrepareGARSOSData.m`:

```matlab
imageSize = [224 224 96];
```

unless an explicit `griddingParams.imageSize` is provided.

The reference geometry is resized with:

```matlab
referenceVolumeC = centeredArrayVolumeResize( ...
    referenceVolumeA, imageSize);
```

at lines 96–97.

`centeredArrayVolumeResize`:

```matlab
B = centeredResize(a.A, newArraySize);
```

and updates the origin:

```matlab
r0 = a.r0 + a.R * diag(a.v) * (cOld - cNew)';
```

at `centeredArrayVolumeResize.m`.

It changes the reference volume’s sampled image grid and origin metadata. It does not transform the 58 acquired kz samples.

## 4.2 Sensitivity-map geometry

The sensitivity maps are first made periodic and resized:

```matlab
periodicMaps = centeredArrayVolumeResize( ...
    channelSensitivityMaps, sensitivityMapArraySize(1:3) .* [3 1 3]);
periodicMaps.A = repmat(channelSensitivityMaps.A, [3 1 3 1]);
```

at `loadPrepareGARSOSData.m`.

They are then spatially resampled onto the reference image geometry:

```matlab
alignedMaps = resampleArrayVolumeToSame( ...
    periodicMaps, referenceVolumeC, 'linear', 0);
```

at line 97.

`resampleArrayVolumeToSame` evaluates source maps at the physical coordinates implied by the reference volume’s:

- rotation `R`;
- voxel size `v`;
- origin `r0`;
- target dimensions.

It does not perform a Fourier transform along the acquired partition axis. It only creates values on the target 224 × 224 × 96 physical grid.

## 4.3 Do map slices and reconstructed slices refer to the same physical plane?

The code is designed for them to do so:

- `alignedMaps.A(:,:,z,:)` is stored on the `referenceVolumeC` 224 × 224 × 96 grid;
- `reconstructedVolume(:,:,z)` is the 224 × 224 × 96 output grid created by `GpuReconstructur20181008`;
- both use the same target `imageSize`;
- the reconstructed volume is centered/resized to that same target grid.

However, the repository does not provide a formal physical metadata object attached to `reconstructedVolume`; it is a numeric array. Exact physical correspondence is therefore an implementation convention, supported by the common image size and centering, rather than something checked at runtime.

The final z indices represent target image-grid positions, not the original 58 partition indices. Their physical coordinates are determined by the target reference geometry:

\[
r(z)
=
r_0 + R\,\operatorname{diag}(v)
\begin{bmatrix}
x\\y\\z
\end{bmatrix}.
\]

The exact numerical z positions require the runtime values of `referenceVolumeA.R`, `referenceVolumeA.v`, and `referenceVolumeA.r0`.

## 4.4 z-axis reversal

No explicit z reversal occurs in `GpuReconstructur20181008` or `shiftAndResize`.

Potential axis orientation changes occur earlier through:

```matlab
frameTrajectory = frameCoordinates * coordinateTransform;
```

in `reconstructGARFrameForSTM.m`.

The geometry transform `CA` comes from the `Cs` construction in `loadPrepareGARSOSData.m`, which includes the `L2` sign-flip matrix:

```matlab
L2 = [diag([1 -1 -1]), ...];
```

Therefore, any physical-axis reversal is attributable to the coordinate transform, not to the final 58-to-96 extraction.

---

# 5. Why the custom 58-point IFFT differs

The custom path does:

```matlab
data4 = reshape(sampledData, ...
    [nReadout, 58, nSpokes, nCoils]);

dataZ4 = fftshift(ifft(ifftshift(data4, 2), [], 2), 2) ...
    * sqrt(58);

selectedZ = dataZ4(:, zIndex, ...);
```

This differs from the working path in several demonstrated ways.

## 5.1 It treats 58 entries as a regular Fourier axis

The custom IFFT assumes:

```text
partition index = regular Fourier-sample index
```

But the working path uses the actual 3-D coordinates:

```matlab
kGrid = sampleCoordinates .* gridSize;
```

and then performs nonuniform kernel deposition.

The 58 entries therefore represent 58 measured trajectory locations, not necessarily 58 equally spaced samples for a direct 58-point DFT.

## 5.2 It uses the wrong transform length

The custom IFFT length is 58.

The working path’s inverse FFT length is the z dimension of `gridSize`:

```text
gridSize(3) = 132
```

not 58.

The final 96 slices are obtained after the 132-point inverse FFT and centered extraction.

## 5.3 It omits 3-D kernel interpolation

The custom path performs no equivalent of:

```matlab
II = ...
DK = ...
accumarray(...)
```

Therefore it does not interpolate each sample onto nearby z-grid bins. It assumes each partition entry already corresponds to one exact regular z-frequency bin.

The working path uses a `[7 7 7]` kernel, including seven-point interpolation in z.

## 5.4 It omits oversampled Cartesian gridding

The working path uses:

```text
z grid: 132
```

from the oversampling factor `1.375`.

The custom path uses only 58.

## 5.5 It selects index 50 from a different coordinate system

In the custom path:

```matlab
zIndex = 50
```

means the 50th sample of a 58-point transformed array.

In the working path:

```matlab
zIndex = 50
```

means the 50th sample of the final centered 96-point image grid after:

```text
nonuniform 3-D gridding
→ 132-point IFFT
→ centered extraction to 96
```

These are different coordinate systems.

## 5.6 The custom path does not replicate `shiftAndResize`

The working path uses:

```matlab
shiftAndResize(..., [224 224 96], 'eee', 'mmm')
```

through the call in `GpuReconstructur20181008`.

The custom path directly indexes the 58-point result and performs no equivalent centered extraction from a 132-point oversampled grid.

## 5.7 It may use a different coordinate transform

The working path transforms sample coordinates before gridding:

```matlab
frameTrajectory = frameCoordinates * coordinateTransform;
```

The custom path’s shown code reshapes the data and applies a partition-axis IFFT, but does not itself perform the same nonuniform coordinate transform and z placement.

Even if the raw 58 entries are ordered correctly, the resulting image-space z index is not guaranteed to correspond to the physical target grid.

---

# 6. Can the existing reconstructor be reused for coil-resolved STM data?

## Approach 1: Directly reuse the existing reconstructor

### Result

```text
224 × 224 × 96
```

### Coil preservation

No. The class sums coils internally.

### In-plane radial sampling

No. The class grids all three dimensions, including x and y. It does not return the original radial readout/spoke organization.

### Physical 96-slice grid

Yes, for the combined image output.

### Assessment

Not suitable directly for STM input.

---

## Approach 2: Modify or parameterize it to return per-coil volumes

### Required change

Instead of accumulating:

```matlab
imageData = imageData + A;
```

store:

```text
imageData(:, iChannel) = A;
```

and reshape to:

```text
[224 224 96 nCoils]
```

The current per-channel loop already performs:

1. 3-D kernel gridding;
2. 3-D `ifftn`;
3. `shiftAndResize`;
4. per-channel apodization.

Therefore this is technically feasible.

### Result

```text
224 × 224 × 96 × nCoils
```

### In-plane radial sampling

No. The result is Cartesian image data. If an in-plane FFT is applied afterward, it produces Cartesian k-space, not the original radial samples.

### Physical 96-slice grid

Yes.

### Mathematical limitation

This output could support image-domain SENSE or coil-resolved Cartesian calibration, but it is not directly the `nReadout × nSpokesPerFrame × nCoils` data expected by the radial STM operator.

### Assessment

Useful if the intended STM formulation is changed to operate on Cartesian coil data, but not a direct solution for the current radial STM reconstruction.

---

## Approach 3: Reuse only the gridding/z-transform internals

### Required operation

Preserve the input arrangement:

```text
[nReadout × 58 × nSpokes × nCoils]
```

and apply only the z-related processing while retaining:

- readout index;
- spoke index;
- coil index;
- in-plane `kx + i ky`.

A z-only equivalent would need to:

1. use each sample’s actual transformed `kz`;
2. grid along a 132-point z Cartesian frequency grid;
3. apply the z-direction inverse FFT;
4. perform the equivalent centered extraction to 96 z positions;
5. select one target z slice.

### Result

For one target z plane:

```text
nReadout × nSpokes × nCoils
```

### In-plane radial sampling

Yes.

### Physical 96-slice grid

Yes, if the z gridding and centered extraction exactly reproduce the z portion of `GpuReconstructur20181008`.

### Important limitation

The current `GpuReconstructur20181008` implementation does not expose a z-only mode. Its `II` and `DK` arrays are built for all three dimensions, and its `grid` method always performs 3-D gridding and `ifftn`.

### Assessment

This is the most mathematically appropriate reuse strategy, but it requires extracting or reimplementing the z-direction part of the gridding pipeline.

---

## Approach 4: Implement an equivalent coil-resolved z operation separately

### Required operation

Implement a dedicated z interpolation/reconstruction path using:

- the same transformed `frameTrajectory(:,3)`;
- the same z grid size `gridSize(3) = 132`;
- the same kernel size `kernelSize(3) = 7`;
- equivalent z kernel weights;
- the same inverse FFT convention;
- the same centered 132-to-96 extraction;
- no coil sum.

### Result

```text
nReadout × nSpokesPerFrame × nCoils
```

for one physical z slice, or:

```text
nReadout × nSpokesPerFrame × nCoils × 96
```

if all target slices are retained.

### In-plane radial sampling

Yes.

### Physical 96-slice grid

Yes, if the z coordinate and centering conventions are matched exactly.

### Assessment

This is the smallest reliable design for the current STM interface, provided the implementation shares the same coordinate transform and z-grid convention as the working path.

---

# 7. Coil-combination boundary and map semantics

The current call in `reconstructGARFrameForSTM.m` is:

```matlab
reconstructor = GpuReconstructur20181008( ...
    gridSize, imageSize, p.Results.KernelSize, ...
    gpuFrameTrajectory, gpuSensitivityArray);
```

The fifth argument is consumed as `coilCombinationProfiles` in `GpuReconstructur20181008.m`.

This means the current working calibration path combines coils before returning `reconstructedVolume`.

With the new forward sensitivity maps, those maps should not be passed directly as combination profiles. The calibration path needs separate combination profiles, such as:

\[
P_c =
\frac{S_c^*}{\sum_j |S_j|^2+\epsilon}.
\]

The radial reconstruction path should continue using the forward maps unchanged:

```matlab
coilImages = x .* senseMaps;
```

in `buildGARSOSRadialEncodingOperator.m`.

This is separate from the 58-to-96 issue, but both issues affect the calibration/reconstruction consistency.

---

# 8. Recommended coil-resolved STM approach

## Recommended design

Do not feed the output of the full `GpuReconstructur20181008` directly into `STM_reconstruction_GARSOS`.

Instead, create a coil-preserving z preprocessing stage equivalent to:

```text
prepared sampledData
    [nReadout × 58 × nSpokes × nCoils]

→ apply the same coordinate transform used by the working reconstructor

→ nonuniform z gridding onto z grid length 132,
   using the same z kernel and z trajectory coordinates

→ inverse FFT along z grid length 132

→ centered extraction/resizing from 132 to 96

→ select target z slice

    [nReadout × nSpokes × nCoils]
```

Then frame the spoke dimension as:

```text
[nReadout × nSpokesPerFrame × nCoils]
```

and retain the corresponding in-plane trajectory:

```text
[nReadout × nSpokesPerFrame]
```

This directly matches the radial STM operator’s expected structure.

## Required consistency points

The z-only path must match:

- `coordinateTransform = prepared.CA`;
- the transformed `kz` values;
- `gridSize(3) = 132`;
- `imageSize(3) = 96`;
- `kernelSize(3) = 7`;
- `ifft` normalization;
- centered shift conventions;
- the 132-to-96 `shiftAndResize` behavior;
- partition/spoke ordering.

It must not combine coils.

---

# 9. Validation test

A useful validation should compare the working combined reconstruction against a coil-preserving reconstruction at several z positions.

Use at least:

```text
z = 1, 24, 50, 72, 96
```

or positions selected from the target physical volume.

## Step A: Run the working reconstructor

For the same frame and same weighted data, obtain:

```matlab
reconstructedVolume(:,:,z)
```

from `GpuReconstructur20181008`.

## Step B: Run the coil-preserving equivalent

Obtain:

```text
coilVolume(:,:,z,c)
```

or, for a radial-preserving z-only path:

```text
coilRadialData(:, :, c, z)
```

## Step C: Apply the same combination profiles

For each z:

\[
x_{\mathrm{combined}}(x,y,z)
=
\sum_c P_c(x,y,z)\,
x_c(x,y,z).
\]

Compare this with:

```matlab
reconstructedVolume(:,:,z)
```

using relative error:

\[
\frac{
\|x_{\mathrm{combined}} - x_{\mathrm{working}}\|_2
}{
\|x_{\mathrm{working}}\|_2
}.
\]

This verifies whether the coil-preserving path reproduces the working reconstruction before using it for STM.

## Step D: Validate the z mapping

For each selected z index, compare:

- the physical coordinate derived from `referenceVolumeC`;
- the corresponding sensitivity-map slice;
- the corresponding combined reconstructed slice.

The same target z index should produce the same physical plane in both arrays.

## Step E: Validate radial preservation

For the proposed STM input, verify:

```text
size(data) = [nReadout, nSpokesPerFrame, nCoils]
```

and compare its in-plane trajectory with:

```matlab
kxy = preparedSlice.kxy;
```

The z-processing step must not alter the `kx`, `ky`, readout, or spoke indexing.

---

# 10. Evidence and unresolved questions

## Directly established by the code

1. `gridSize` is `ceil(1.375 .* imageSize)`, giving `[308 308 132]` for `[224 224 96]`.
2. The reconstructor deposits samples onto the 3-D Cartesian grid using kernel neighborhoods.
3. The gridded array is inverse transformed with `ifftn`.
4. `shiftAndResize` converts the 308 × 308 × 132 result to 224 × 224 × 96.
5. The 58 input partition locations are not directly transformed by a 58-point IFFT in the working path.
6. Coil images are gridded separately but summed inside `GpuReconstructur20181008`.
7. The current reconstructor cannot return coil-resolved output without modification.
8. `prepareGARSOSSliceDataForMCNUFFT` performs a direct 58-point partition IFFT and restricts `zIndex` to `1:nPartition`; this is a different path from the working 3-D nonuniform gridding path.

## Not established without runtime inspection

1. The exact numerical values, spacing, and range of `unique(frameTrajectory(:,3))`.
2. Whether the transformed z coordinates are uniformly spaced after multiplication by `CA`.
3. The exact physical coordinates of output slices 1 through 96, because runtime `R`, `v`, and `r0` values are required.
4. Whether the input partition zero-padding and trajectory construction are perfectly consistent for the specific acquisition.
5. The numerical error between a proposed z-only coil-preserving implementation and the existing 3-D reconstructor.

## Most important unresolved implementation issue

The current custom code assumes that the 58 partition samples are a regular Fourier axis and that index 50 of a 58-point IFFT is target image slice 50. The working code demonstrates that the target image slices instead arise from a 3-D nonuniform gridding and oversampled inverse-transform pipeline.

Therefore, `zIndex = 50` in the custom 58-point path is not equivalent to `zIndex = 50` in the working 96-slice reconstruction.