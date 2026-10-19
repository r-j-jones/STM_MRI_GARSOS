# PISCO 3D STM prototype: Phases C, D, and preliminary E

## New files

- `nullspace_vectors_G_matrix_3D_STM_v2.m`: SVD or block orthogonal iteration, optional Rayleigh–Ritz, optional 3D FFT interpolation, phase alignment, eigenvalue interpolation, and map normalization.
- `PISCO_sensitivity_map_estimation_3D_STM.m`: PISCO-style 3D wrapper preserving C-nullspace and G construction while routing the final eigensolver to the STM function.
- `example_3D_STM.m`: modified PISCO example illustrating SVD and orthogonal-iteration calls.
- `test_nullspace_vectors_G_matrix_3D_STM_v2.m`: synthetic Phase A-D tests.

## Phase C: Rayleigh–Ritz

The iterative solver exposes:

```matlab
'RayleighRitz',1
```

After block orthogonal iteration, it forms the projected matrix `H=S'*A*S`, diagonalizes each `L x L` page, sorts its eigenvalues ascending, and rotates `S` into approximate individual eigenvectors. With the option disabled, it returns the orthonormal subspace basis and diagonal Rayleigh quotients, which need not be individually ordered eigenvectors.

## Phase D: FFT interpolation

The implementation interpolates only spatial dimensions 1-3 and preserves dimensions 4 (`Nc`) and 5 (`L`). Options are:

- `'FFT_interpolation',1|0`
- `'gauss_win_param',100`
- `'PhaseAlignment','calibration'|'first-channel'|'none'`
- `'InterpolateEigenvalues',1|0`
- `'NormalizeMaps',1|0`

Calibration phase alignment is performed independently for each STM index `L`. It fixes complex scalar phase but cannot remove arbitrary unitary rotations within clustered eigenspaces.

## Preliminary Phase E

`PISCO_sensitivity_map_estimation_3D_STM.m` preserves the PISCO stages:

1. 3D C-nullspace calculation.
2. 3D G-matrix construction.
3. Multi-vector 3D STM eigensolver.

It assumes the fourth dimension of `kCal` is the virtual-channel/time dimension. It is a prototype and must not be used on arrays containing both receive-coil and time dimensions until that data contract is explicitly resolved.

## Required path setup

The folder containing the new files and the PISCO repository must be on the MATLAB path. The PISCO package folder `+utils` must resolve correctly. Avoid recursively adding multiple repositories that contain duplicate function names.

## Suggested order of use

1. Run `test_nullspace_vectors_G_matrix_3D_STM_v2`.
2. Run `example_3D_STM` with SVD on a small calibration block.
3. Compare orthogonal iteration to SVD using subspace projectors.
4. Verify the GARSOS calibration dimension order before substituting real data.
5. Integrate into the GARSOS reconstruction only after validating output dimensions, residuals, memory use, and downstream reconstruction.