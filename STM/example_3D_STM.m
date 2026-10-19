% EXAMPLE_3D_STM
% Prototype adaptation of PISCO example_3D.m for L three-dimensional STMs.
% It assumes kData is arranged [kx,ky,kz,time/virtual-channel].
% Verify this interpretation before using physical receive-coil data.

clear; close all; clc
load('./data/3D_GRE_phantom_data.mat') % Replace with GARSOS 3D kCal source.

[N1,N2,N3,Nvirtual]=size(kData);
L=2;
if L>Nvirtual, error('L must not exceed the virtual-channel/time dimension.'); end

% Central 3D calibration block.
cal_length=32;
cx=ceil(N1/2)+utils.even_pisco(N1);
cy=ceil(N2/2)+utils.even_pisco(N2);
cz=ceil(N3/2)+utils.even_pisco(N3);
ix=cx+(-floor(cal_length/2):floor(cal_length/2)-utils.even_pisco(cal_length/2));
iy=cy+(-floor(cal_length/2):floor(cal_length/2)-utils.even_pisco(cal_length/2));
iz=cz+(-floor(cal_length/2):floor(cal_length/2)-utils.even_pisco(cal_length/2));
kCal=kData(ix,iy,iz,:);
dim_stm=[N1 N2 N3];

% Shared PISCO options.
common={ ...
 'tau',3,'threshold',0.08,'kernel_shape',1, ...
 'FFT_nullspace_C_calculation',1,'FFT_interpolation',1,'interp_zp',24, ...
 'gauss_win_param',100,'sketched_SVD',1,'sketch_dim',500, ...
 'visualize_C_matrix_sv',0,'PhaseAlignment','calibration', ...
 'InterpolateEigenvalues',1,'NormalizeMaps',0,'verbose',1};

% Phase A reference: paged SVD plus Phase D interpolation.
[ST_svd,eig_svd,diag_svd]=PISCO_sensitivity_map_estimation_3D_STM( ...
 kCal,dim_stm,L,common{:},'Solver','svd','RayleighRitz',1);

% Phases B/C/D: block orthogonal iteration, Rayleigh-Ritz, interpolation.
[ST_oi,eig_oi,diag_oi]=PISCO_sensitivity_map_estimation_3D_STM( ...
 kCal,dim_stm,L,common{:},'Solver','orthogonal-iteration', ...
 'M',100,'Tolerance',1e-8,'CheckEvery',5,'RandomSeed',0,'RayleighRitz',1);

fprintf('SVD output: %s\n',mat2str(size(ST_svd)));
fprintf('OI output:  %s\n',mat2str(size(ST_oi)));
fprintf('OI iterations: %d; converged: %d\n',diag_oi.iterations,diag_oi.converged);

% Compare subspaces at one voxel. The Nc-by-L columns may differ by phase
% or unitary rotation, so compare projectors rather than raw vectors.
voxel=round(dim_stm/2);
Vs=reshape(ST_svd(voxel(1),voxel(2),voxel(3),:,:),[Nvirtual L]);
Vo=reshape(ST_oi(voxel(1),voxel(2),voxel(3),:,:),[Nvirtual L]);
projectorDifference=norm(Vs*Vs'-Vo*Vo','fro');
fprintf('Center-voxel projector difference: %.3e\n',projectorDifference);

% Display magnitudes for one z slice. Each row/column montage arrangement
% depends on utils.mdisp from PISCO.
slc=voxel(3);
for ell=1:L
 figure;
 imagesc(utils.mdisp(abs(squeeze(ST_oi(:,:,slc,:,ell)))));
 axis image tight off; colormap gray;
 title(sprintf('3D STM magnitude, L index %d, z slice %d',ell,slc));
end

figure;
imagesc(utils.mdisp(abs(squeeze(eig_oi(:,:,slc,:)))));
axis image tight off; colormap gray; colorbar;
title('Selected 3D STM eigenvalue maps');

% Optional alternative experiments:
%   'FFT_interpolation',0 requires the requested dimensions to match the G grid.
%   'RayleighRitz',0 returns an orthonormal subspace basis but not guaranteed
%       individually ordered Ritz eigenvectors.
%   'PhaseAlignment','first-channel' or 'none'.
%   'NormalizeMaps',1 normalizes each L map across virtual channels.