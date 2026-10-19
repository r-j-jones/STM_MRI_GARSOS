# Adapting PISCO 3D to compute multiple 3D STMs

## Purpose

This document outlines how to adapt PISCO's three-dimensional local `G`-matrix eigensolver to compute `L` three-dimensional spatiotemporal maps (STMs). It is based on comparing:

- STM-2D `nullspace_vectors_G_matrix_2D.m`, which returns `L` STMs;
- PISCO-2D `nullspace_vectors_G_matrix_2D.m`, which returns one sensitivity-map vector;
- PISCO-3D `nullspace_vectors_G_matrix_3D.m`, which returns one sensitivity-map vector.

The proposed implementation combines PISCO-3D's spatial organization with STM-2D's multi-vector eigensolver. In the STM interpretation, PISCO's channel dimension is treated as a virtual-channel/time-frame dimension rather than a physical receive-coil dimension.

## Main conclusion

The adaptation is feasible and localized. At each 3D spatial location, PISCO already constructs an `Nc x Nc` Hermitian positive-semidefinite matrix `G(x,y,z)`. PISCO currently estimates one vector associated with its smallest eigenvalue. A 3D STM implementation must instead estimate the `L` vectors associated with the `L` smallest eigenvalues.

The natural output shape is:

```matlab
[N1_g, N2_g, N3_g, Nc, L]
```

where:

- dimensions 1-3 are spatial;
- dimension 4 contains the `Nc` virtual-channel/time-frame coefficients;
- dimension 5 indexes the `L` STMs.

The `Nc` and `L` dimensions should not be flattened together.

## Input and output conventions

### Input calibration data

For direct 3D STM estimation:

```matlab
kCal: N1_cal x N2_cal x N3_cal x Nc
```

For the intended GARSOS use, `Nc` represents temporal frames or virtual channels.

### Input local matrices

```matlab
G: N1_g x N2_g x N3_g x Nc x Nc
```

At voxel `(i,j,k)`, `G(i,j,k,:,:)` is the local `Nc x Nc` matrix.

### Output maps

Before spatial interpolation:

```matlab
ST_maps: N1_g x N2_g x N3_g x Nc x L
```

A future interpolation phase may produce:

```matlab
ST_maps: N1 x N2 x N3 x Nc x L
```

### Output eigenvalues

The implementation in this folder returns only the selected eigenvalues:

```matlab
eigenVal: N1_g x N2_g x N3_g x L
```

They are ordered from smallest to largest at each voxel.

## Relationship between the algorithms

### Direct SVD

Both existing 2D and 3D functions reshape spatially indexed `G` matrices into paged matrices:

```matlab
Nc x Nc x numberOfVoxels
```

MATLAB `pagesvd` returns singular values in descending order. Because `G` is expected to be Hermitian positive semidefinite, its singular vectors are its eigenvectors, and the final columns correspond to its smallest eigenvalues.

PISCO-3D selects one vector:

```matlab
Vpage(:,end,:)
```

STM-2D selects `L` vectors:

```matlab
Vpage(:,Nc-L+1:Nc,:)
```

For a consistent STM interface, the selected vectors and values should then be reordered so the returned STM index corresponds to ascending eigenvalue.

### Power iteration versus orthogonal iteration

PISCO forms:

```matlab
G_null = I - G/patchSize
```

The smallest eigenvalue of normalized `G` maps to the largest eigenvalue of `G_null`. Ordinary power iteration therefore estimates one desired vector and is appropriate for `L=1`.

For `L>1`, STM-2D uses a block method:

1. Initialize `L` vectors per voxel.
2. Multiply the whole block by `G_null`.
3. Orthonormalize the `L` columns.
4. Repeat for `M` iterations.

This is block power iteration, also called orthogonal or subspace iteration. It estimates the dominant `L`-dimensional invariant subspace of `G_null`, which is the subspace associated with the `L` smallest eigenvalues of `G`.

PISCO's second vector `S2` is used only as a convergence diagnostic for its single-vector power iteration. It is not a general `L`-vector solver and should not be reused as though it were one.

## Phase A: exact paged-SVD reference

The first implementation should use the direct SVD path with no FFT interpolation.

Steps:

1. Validate `1 <= L <= Nc`.
2. Reshape `G` from:

   ```matlab
   [N1_g,N2_g,N3_g,Nc,Nc]
   ```

   to:

   ```matlab
   [Nc,Nc,N1_g*N2_g*N3_g]
   ```

3. Hermitian-symmetrize each page to suppress small numerical asymmetry:

   ```matlab
   G = 0.5*(G + G');
   ```

4. Run paged SVD.
5. Select the final `L` right singular vectors.
6. Reverse their within-selection order so output eigenvalues are ascending.
7. Normalize eigenvalues by `patchSize`.
8. Return maps with shape `[N1_g,N2_g,N3_g,Nc,L]`.

This provides a reference solution for validating the iterative solver.

## Phase B: block orthogonal iteration

The iterative method starts from normalized:

```matlab
A = G/patchSize
```

and constructs a shifted matrix whose dominant subspace matches the smallest-eigenvalue subspace of `A`.

The original STM and PISCO code uses `I-A`. The implementation supplied here uses a safe per-page shift:

```matlab
B = alpha*I - A
```

where `alpha` is at least the largest eigenvalue of each page. This retains the desired eigenvectors and avoids relying on the assumption that all eigenvalues of `A` are bounded by one.

Steps:

1. Initialize `S` with shape:

   ```matlab
   [Nc,L,numVoxels]
   ```

2. Orthonormalize each page.
3. Repeat `M` times:

   ```matlab
   S = B*S
   S = pagewise_orthonormalize(S)
   ```

4. Perform a Rayleigh-Ritz refinement:

   ```matlab
   H = S'*A*S
   ```

5. Diagonalize each small `L x L` matrix `H`.
6. Rotate `S` by those projected eigenvectors.
7. Sort the selected eigenvalues in ascending order.
8. Return the rotated, ordered vectors as the `L` STMs.

### Why Rayleigh-Ritz is useful

Block orthogonal iteration primarily converges to a subspace. Its columns need not themselves be individually ordered eigenvectors. Rayleigh-Ritz rotates the converged basis within the subspace so that each returned column is an approximate eigenvector with a corresponding ordered eigenvalue.

This matters when comparing residuals or interpreting each `L` map separately.

## Convergence considerations

PISCO's single-vector convergence ratio is based on the separation of the first and second dominant modes of the shifted matrix. For an `L`-dimensional subspace, convergence depends on separation between the `L`th and `(L+1)`th modes.

The supplied Phase B implementation therefore uses a fixed iteration count and reports a subspace residual when requested:

```matlab
R = A*S - S*(S'*A*S)
```

For validation, compare the resulting subspace to the SVD reference using projectors:

```matlab
norm(S_svd*S_svd' - S_oi*S_oi','fro')
```

Projector comparison is preferable to direct vector comparison because eigenvectors have arbitrary complex phase, and a repeated or tightly clustered eigenvalue allows arbitrary unitary rotations within its invariant subspace.

## FFT interpolation: intentionally deferred

The supplied implementation covers Phase A and Phase B only. It requires:

```matlab
FFT_interpolation = 0
```

PISCO-3D's interpolation assumes maps shaped as:

```matlab
[N1_g,N2_g,N3_g,Nc]
```

A 3D STM implementation must explicitly preserve the extra `L` dimension:

```matlab
[N1_g,N2_g,N3_g,Nc,L]
```

The future interpolation window should cover both trailing dimensions while transforming only dimensions 1-3. Phase alignment would need to be computed separately for each `L` map. Independent phase alignment does not resolve arbitrary unitary rotations when selected eigenvalues are repeated or nearly repeated, so spatial subspace alignment may eventually be necessary.

## Implementation files

### `nullspace_vectors_G_matrix_3D_STM.m`

Implements:

- Phase A paged-SVD reference;
- Phase B block orthogonal iteration;
- pagewise modified Gram-Schmidt reorthogonalization;
- Rayleigh-Ritz extraction;
- ascending eigenvalue order;
- optional deterministic random seed;
- optional residual reporting;
- no FFT interpolation.

The default solver is `svd` so the first use favors a reliable reference. Select the iterative method with:

```matlab
'Solver','orthogonal-iteration'
```

### `test_nullspace_vectors_G_matrix_3D_STM.m`

Constructs synthetic Hermitian positive-semidefinite `G` matrices with known spectra and verifies:

- output dimensions;
- ascending eigenvalue order;
- orthonormality;
- eigen residuals;
- agreement of the iterative and SVD subspace projectors;
- agreement with MATLAB `eig` at one voxel;
- `L=1` dimensional behavior.

## Suggested calling pattern

### Phase A

```matlab
[ST_svd,eig_svd] = nullspace_vectors_G_matrix_3D_STM( ...
    kCal,N1,N2,N3,L,G,patchSize, ...
    'Solver','svd', ...
    'FFT_interpolation',0);
```

### Phase B

```matlab
[ST_oi,eig_oi] = nullspace_vectors_G_matrix_3D_STM( ...
    kCal,N1,N2,N3,L,G,patchSize, ...
    'Solver','orthogonal-iteration', ...
    'M',100, ...
    'Tolerance',1e-8, ...
    'FFT_interpolation',0, ...
    'verbose',1);
```

## Validation priorities before GARSOS integration

1. Run the supplied synthetic test.
2. Test `L=1` against PISCO-3D up to complex phase.
3. Test `L>1` against the SVD reference using subspace projectors.
4. Confirm GARSOS `kCal` is ordered as `[kx,ky,kz,time]` before treating time as the virtual-channel dimension.
5. If physical receive coils and time are both present, do not silently flatten them without defining the intended joint-channel model.
6. Test memory use on realistic dimensions. Paged SVD and full materialization of `G` can be expensive.
7. Add 3D FFT interpolation only after low-resolution eigensolver results are correct.

## Feasibility assessment

Phase A is a small extension of PISCO-3D: add `L`, select the final `L` singular vectors, reorder them, and reshape with a fifth dimension.

Phase B is also practical: transfer STM-2D's block orthogonal iteration to the 3D page organization and add a Rayleigh-Ritz refinement. The spatial dimensionality changes the number of pages and output reshape, but not the local `Nc x Nc` eigensolver mathematics.

The areas requiring the most care beyond these phases are convergence diagnostics, memory consumption, interpolation over the extra `L` dimension, and spatial alignment of bases when eigenvalues are clustered.