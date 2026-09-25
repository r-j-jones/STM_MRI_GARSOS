# GAR-SOS Pipeline Refactor Report

## Agent Activity Log

### Initial Plan

> I’ll implement the shared loader and refactor the three pipeline entry points together, then audit the cached-manifest validation, usage text, and `_v1` references before reporting the workflow and call tree.

### Initial Searches and Environment Check

- Searched for `reconstructGARFrameForSTM_v1`.
- Searched for files matching `GARSOS_MotionModeling/*Loader*.m`.
- Ran:

  ```bash
  command -v matlab || true
  command -v octave || true
  ```

### Agent Read Failure

- **Action:** Read Agent
- **Agent ID:** `5011ada4-8415-4ab-be36-8028aca67fe8`
- **Wait:** `true`
- **Timeout:** 180 seconds
- **Result:** Failed

### Source Files Reviewed

- Read lines 1–120 of:
  `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/prepareGARSOSCalibrationData.m`
- Read lines 1–140 of:
  `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/prepareGARSOSSTMData.m`
- Read lines 1–130 of:
  `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/main_STM_GARSOS.m`

### Refactor Decision

> The requested refactor is feasible. I’m proceeding with a single-sequence shared preparation struct and an optional prepared-data handoff so calibration and STM stages reuse the same raw/preprocessed arrays instead of loading them twice.

### Files Edited

- `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/loadPrepareGARSOSData.m`
- `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/prepareGARSOSSTMData.m`
- `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/prepareGARSOSCalibrationData.m`
- `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/main_STM_GARSOS.m`
- `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/reconstructGARFrameForSTMCoordinateAware.m`
- `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/radial_vibe_18_5_7_7_imFIAT_permuted.m`
- `/Users/robertjones/Desktop/hmm_code_matlab/GARSOS_MotionModeling/runGARSOSSTMIntegration.m`

### Failed Patch Attempt

An initial patch to `prepareGARSOSCalibrationData.m` failed. The intended changes included:

- Changing the function return values from:

  ```matlab
  [manifest, manifestFile, calibData, stmCalibOutputs]
  ```

  to:

  ```matlab
  [manifest, manifestFile, stmCalibOutputs]
  ```

- Adding optional name/value handling for:

  ```matlab
  'PreparedData', preparedData
  ```

- Validating cached sequence-level `calibData`.
- Replacing duplicated raw-data preparation with:

  ```matlab
  loadPrepareGARSOSData(...)
  ```

- Reusing fields from the returned `preparedData` struct.

After the patch failed, the file was edited through subsequent operations.

### Additional Reviews, Searches, and Checks

The following source ranges were reviewed:

- Lines 110–225 of `prepareGARSOSCalibrationData.m`
- Lines 395–455 of `prepareGARSOSCalibrationData.m`
- Lines 115–265 of `prepareGARSOSCalibrationData.m`
- Lines 1–120 of `prepareGARSOSSTMData.m`

Searches were performed for:

- `calibData`
- `function [`
- `reconstructGARFrameForSTM_v1`
- `localDisplayUsage`
- `Reruun`
- `RerunWithoutSaving`
- `radial_vibe_18_5_7_7_imFIAT_permuted`
- `channelSensitivityMaps`
- `r0Machine`
- `preparedData`

The environment and recently modified files were also checked with commands such as:

```bash
command -v matlab || true
command -v octave || true
find GARSOS_MotionModeling -maxdepth 1 -type f -newermt '45 minutes ago' -print | sort
```

---

## Implementation Summary

The shared GAR-SOS preparation refactor was implemented, and the active STM pipeline was updated.

### Files Added

- `loadPrepareGARSOSData.m`

### Files Modified

- `main_STM_GARSOS.m`
- `prepareGARSOSSTMData.m`
- `prepareGARSOSCalibrationData.m`
- `runGARSOSSTMIntegration.m`
- `reconstructGARFrameForSTMCoordinateAware.m`
- `radial_vibe_18_5_7_7_imFIAT_permuted.m`

Active references to `reconstructGARFrameForSTM_v1` were changed to `reconstructGARFrameForSTM`.

---

### Shared Loader

`loadPrepareGARSOSData.m` now centralizes:

- Raw Siemens data loading
- `mapVBVD`
- Coil-sensitivity estimation
- Radial preprocessing through `prepareRadialVibeDataForGridding2_3`
- Dimension and image/grid sizing
- `Cs`, `CA`, and `Cb` derivation from `L1`, `L2`, and `L3`
- Periodic sensitivity-map extension
- Scanner-coordinate origin adjustment
- Sensitivity-map resampling to the target image geometry

It returns an organized `prepared` struct containing:

```matlab
prepared.sampledData
prepared.kSpaceLocations
prepared.densityCompensation
prepared.spokeIndex
prepared.timeStamps

prepared.nPartitions
prepared.nSpokes
prepared.nSamples
prepared.nChannels

prepared.sensitivityMaps
prepared.alignedSensitivityMaps
prepared.channelSensitivityMaps

prepared.Cs
prepared.CA
prepared.Cb
prepared.L1
prepared.L2
prepared.L3

prepared.dimensions
prepared.metadata
prepared.rawDataFile
```

---

### `prepareGARSOSCalibrationData`

`prepareGARSOSCalibrationData.m` now:

- Accepts the requested argument order:

  ```matlab
  prepareGARSOSCalibrationData( ...
      metaManifest, reconstructionDir, referenceVolumeA, ...
      frameParams, griddingParams, forceUpdate)
  ```

- Accepts optional name/value input:

  ```matlab
  'PreparedData', prepared
  ```

- Reuses the shared prepared data when supplied.
- Calls `loadPrepareGARSOSData` itself when prepared data are not supplied.
- Preserves the correctly spelled option:

  ```matlab
  'RerunWithoutSaving'
  ```

- Validates cached `manifest.sequences(1).calibData`.
- Throws explicit errors if cached `calibData` is missing or lacks:

  ```matlab
  kCal
  nbSpokesPerFrame
  zSlice
  calibrationSize
  ```

- Retains the `calibData` field in the saved manifest.
- Returns only:

  ```matlab
  [manifest, manifestFile, stmCalibOutputs]
  ```

- Keeps large reconstruction arrays in `stmCalibOutputs`, rather than in the manifest.

The cached path now treats the manifest’s sequence-level `calibData` as authoritative.

---

### `prepareGARSOSSTMData`

`prepareGARSOSSTMData.m` now:

- Accepts the following options:
  - `FrameParams`
  - `GriddingParams`
  - `PreparedData`
- Uses `loadPrepareGARSOSData` if prepared data are not supplied.
- Calls `prepareGARSOSCalibrationData`.
- Passes the shared prepared data into calibration.
- Stores the calibration manifest and calibration data in the returned struct.
- Avoids calling `radial_vibe_18_5_7_7_imFIAT_permuted`.

The returned struct now includes:

```matlab
prepared.calibrationManifest
prepared.calibrationManifestFile
prepared.calibData
prepared.stmCalibOutputs
prepared.frameParams
prepared.griddingParams
```

---

### `main_STM_GARSOS`

`main_STM_GARSOS.m` now:

1. Loads the Fessler nonuniform fast Fourier transform toolbox.
2. Validates input paths and options.
3. Loads the reference DICOM volume.
4. Runs `radial_vibe_18_5_4_timestamp_imFIAT`.
5. Constructs:

   ```matlab
   frameParams
   griddingParams
   ```

6. Calls `loadPrepareGARSOSData` once.
7. Passes the resulting prepared data into `prepareGARSOSSTMData`.
8. Uses the prepared calibration data when constructing the STM integration.
9. Runs `runGARSOSSTMIntegration`.
10. Returns the final result with prepared data, calibration data, the radial operator, and the STM operator.

The default frame configuration remains configurable through:

```matlab
opts.nbSpokesPerFrame
opts.zSlice
opts.calibrationSize
```

---

### Avoiding Duplicate Calibration Reconstruction

`runGARSOSSTMIntegration.m` now checks for:

```matlab
prepared.calibData.kCal
```

When that field is present, the function reuses it instead of regenerating `kCal` through another conventional reconstruction.

---

## Main Workflow

```text
main_STM_GARSOS
│
├── Add Fessler NUFFT toolbox
│
├── Validate input arguments and paths
│
├── Load reconstruction/input.mat
│   ├── load kSpaceID1
│   └── load kSpaceID2
│
├── loadReferenceDicomImage
│
├── parseKspaceIds
│
├── radial_vibe_18_5_4_timestamp_imFIAT
│   └── generate/load acquisition metadata manifest
│
├── Construct frameParams
│   ├── nSpokesPerFrame
│   ├── zSlice
│   └── calibrationSize
│
├── Construct griddingParams
│   ├── UseGpu
│   ├── DisplaySlice
│   ├── ReturnDiagnostics
│   └── gridOversamplingFactor
│
├── loadPrepareGARSOSData
│   ├── readVD11MultiRaidFileStructure
│   ├── mapVBVD
│   ├── estimateCoilSensitivitieMaps2
│   ├── prepareRadialVibeDataForGridding2_3
│   ├── centeredArrayVolumeResize
│   ├── derive L1/L2/L3
│   ├── derive Cs/CA/Cb
│   └── resampleArrayVolumeToSame
│
├── prepareGARSOSSTMData
│   └── prepareGARSOSCalibrationData
│       ├── validate cached manifest/calibData
│       ├── reconstructGARFrameForSTM
│       ├── extract kCal
│       ├── collect diagnostics
│       └── return stmCalibOutputs
│
└── runGARSOSSTMIntegration
    ├── use existing kCal when available
    ├── prepare transformed radial slice data
    ├── validate trajectory/data correspondence
    ├── buildGARSOSRadialEncodingOperator
    │   └── buildGARSOSRadialNUFFTPlans
    │       ├── nufft_init
    │       ├── nufft
    │       └── nufft_adj
    ├── buildGARSOSSTMOperators
    └── return STM/radial reconstruction result
```

---

## Function Call Tree

```text
main_STM_GARSOS
├── loadReferenceDicomImage
├── parseKspaceIds
├── radial_vibe_18_5_4_timestamp_imFIAT
├── loadPrepareGARSOSData
│   ├── readVD11MultiRaidFileStructure
│   ├── mapVBVD
│   ├── estimateCoilSensitivitieMaps2
│   ├── prepareRadialVibeDataForGridding2_3
│   ├── centeredArrayVolumeResize
│   ├── resampleArrayVolumeToSame
│   └── derive Cs / CA / Cb
├── prepareGARSOSSTMData
│   └── prepareGARSOSCalibrationData
│       ├── loadPrepareGARSOSData
│       │   └── skipped when PreparedData is supplied
│       └── reconstructGARFrameForSTM
└── runGARSOSSTMIntegration
    ├── prepareGARSOSSliceDataForMCNUFFT
    ├── validateGARSOSTrajectoryData
    ├── buildGARSOSRadialEncodingOperator
    │   └── buildGARSOSRadialNUFFTPlans
    └── buildGARSOSSTMOperators
```

---

## Remaining Concerns

1. **MATLAB was unavailable.** The end-to-end pipeline could not be run or syntax-checked through MATLAB or Octave.

2. **The implementation preserves the requested single-sequence assumption.** The calibration function still loops over manifest sequences, but the shared handoff and cached path currently use sequence 1.

3. **The main function contains a machine-specific Fessler setup path:**

   ```matlab
   /RadOnc-MRI1/Student_Folder/rjones/toolboxes/fessler/irt/setup.m
   ```

   This path will need to be changed on systems where it does not exist.

4. **Older archival files still contain `_v1` references.** This includes files under `all_HMM`. The active `GARSOS_MotionModeling` pipeline references the renamed `reconstructGARFrameForSTM`, but archival copies were not changed.

5. **The shared loader depends on existing toolboxes.** These include the Siemens, ArrayVolume, GAR-SOS, and nonuniform fast Fourier transform toolboxes.

6. **Runtime validation remains pending.** No MATLAB runtime tests or measured operator-adjointness results could be produced in the available environment.

---

## Delegated Task Audit

### Checkpoint Review

Read:

`/Users/robertjones/.copilot/session-state/d4de9f3d-5dbb-4006-8422-e79882c610b3/checkpoints/001-debugging-gar-sos-reconstructi.md`

The following files were also located and reviewed:

- `loadPrepareGARSOSData.m`
- `prepareGARSOSCalibrationData.m`
- `prepareGARSOSSTMData.m`
- `main_STM_GARSOS.m`
- `reconstructGARFrameForSTM.m`
- `runGARSOSSTMIntegration.m`
- `reconstructGARFrameForSTMCoordinateAware.m`
- `orig__prepareGARSOSCalibrationData.m`

Specific source ranges reviewed included:

- Lines 1–240 of `prepareGARSOSCalibrationData.m`
- Lines 241–520 of `prepareGARSOSCalibrationData.m`
- Lines 120–260 of `prepareGARSOSCalibrationData.m`
- Lines 390–470 of `prepareGARSOSCalibrationData.m`
- Lines 1–260 of `reconstructGARFrameForSTM.m`
- Lines 1–220 of `runGARSOSSTMIntegration.m`
- Lines 1–180 of `reconstructGARFrameForSTMCoordinateAware.m`
- Lines 1–60 of `orig__prepareGARSOSCalibrationData.m`

### Static Checks

Searches were performed for:

```text
radial_vibe_18_5_7_7_imFIAT_permuted
reconstructGARFrameForSTM_v1
_v1
calibData
prepareGARSOSCalibrationData
```

Static shell checks included:

```bash
command -v octave || true
command -v matlab || true
git diff --no-index /dev/null GARSOS_MotionModeling/loadPrepareGARSOSData.m >/dev/null
git diff -- GARSOS_MotionModeling
```

File statistics were collected with:

```bash
python - <<'PY'
from pathlib import Path

root = Path('GARSOS_MotionModeling')

for name in [
    'loadPrepareGARSOSData.m',
    'prepareGARSOSCalibrationData.m',
    'prepareGARSOSSTMData.m',
    'main_STM_GARSOS.m'
]:
    path = root / name
    text = path.read_text()
    function_count = sum(
        1 for line in text.splitlines()
        if line.lstrip().startswith('function ')
    )
    print(
        name,
        'lines',
        text.count('\n') + 1,
        'functions',
        function_count
    )
PY
```

Targeted searches included:

```bash
grep -RInE \
'radial_vibe_18_5_7_7_imFIAT_permuted|reconstructGARFrameForSTM_v1|_v1' \
GARSOS_MotionModeling/{prepareGARSOSCalibrationData.m,prepareGARSOSSTMData.m,main_STM_GARSOS.m,loadPrepareGARSOSData.m,runGARSOSSTMIntegration.m} \
|| true
```

and:

```bash
grep -RIn \
'prepareGARSOSCalibrationData' \
GARSOS_MotionModeling \
--include='*.m'
```

---

## Delegated Agent Result

The background agent `Refactor GAR-SOS pipeline` completed.

- **Agent ID:** `5011ada4-8415-4ab1-be36-8028aca67fe8`
- **Wait:** `false`
- **Result:** Complete

The delegated refactor agent confirmed that the shared GAR-SOS preparation changes were implemented and statically audited. No additional changes were required based on that report.

MATLAB and Octave remained unavailable in the environment, so runtime validation is still pending.
