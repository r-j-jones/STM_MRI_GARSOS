function [ST_maps,eigenValues,diagnostics] = PISCO_sensitivity_map_estimation_3D_STM(kCal,dim_stm,L,varargin)
%PISCO_SENSITIVITY_MAP_ESTIMATION_3D_STM Prototype PISCO-to-STM wrapper.
%
% This preserves PISCO's 3D C-nullspace and G-matrix construction, then
% calls nullspace_vectors_G_matrix_3D_STM_v2 to return L 3D STMs.
%
% Required PISCO utilities on the path:
%   utils.nullspace_vectors_C_matrix_3D
%   utils.G_matrices_3D
%
% Input
%   kCal     N1_cal x N2_cal x N3_cal x Nt. Dimension 4 is interpreted as
%            virtual channels/time frames.
%   dim_stm  [N1 N2 N3] desired spatial dimensions.
%   L        number of STM vectors.
%
% Major options mirror PISCO where applicable:
%   tau, threshold, kernel_shape, FFT_nullspace_C_calculation,
%   FFT_interpolation, interp_zp, gauss_win_param, sketched_SVD,
%   sketch_dim, visualize_C_matrix_sv, verbose.
%
% STM eigensolver options:
%   Solver ('svd' or 'orthogonal-iteration'), M, Tolerance, CheckEvery,
%   RandomSeed, RayleighRitz, PhaseAlignment, InterpolateEigenvalues,
%   NormalizeMaps.

p=inputParser;
p.addRequired('kCal',@(x)isnumeric(x)&&ndims(x)==4);
p.addRequired('dim_stm',@(x)isnumeric(x)&&isvector(x)&&numel(x)==3);
p.addRequired('L',@(x)isnumeric(x)&&isscalar(x)&&x>=1&&x==floor(x));
p.addParameter('tau',3,@isscalar);
p.addParameter('threshold',0.05,@isscalar);
p.addParameter('kernel_shape',1,@isscalar);
p.addParameter('FFT_nullspace_C_calculation',1,@isscalar);
p.addParameter('FFT_interpolation',1,@isscalar);
p.addParameter('interp_zp',24,@isscalar);
p.addParameter('gauss_win_param',100,@isscalar);
p.addParameter('sketched_SVD',1,@isscalar);
p.addParameter('sketch_dim',500,@isscalar);
p.addParameter('visualize_C_matrix_sv',0,@isscalar);
p.addParameter('Solver','orthogonal-iteration');
p.addParameter('M',100,@isscalar);
p.addParameter('Tolerance',1e-8,@isscalar);
p.addParameter('CheckEvery',5,@isscalar);
p.addParameter('RandomSeed',0,@isscalar);
p.addParameter('RayleighRitz',1,@isscalar);
p.addParameter('PhaseAlignment','calibration');
p.addParameter('InterpolateEigenvalues',1,@isscalar);
p.addParameter('NormalizeMaps',0,@isscalar);
p.addParameter('verbose',1,@isscalar);
p.parse(kCal,dim_stm,L,varargin{:});
r=p.Results;

if L>size(kCal,4), error('L must not exceed size(kCal,4).'); end
if r.FFT_interpolation
    gridSize=size(kCal,1:3)+r.interp_zp;
else
    gridSize=dim_stm;
end

if r.verbose
 fprintf('3D STM prototype: constructing C nullspace and G matrices.\n');
 fprintf('Virtual-channel/time dimension: %d; requested L: %d.\n',size(kCal,4),L);
end

totalTimer=tic;
U=utils.nullspace_vectors_C_matrix_3D(kCal, ...
    'tau',r.tau,'threshold',r.threshold,'kernel_shape',r.kernel_shape, ...
    'FFT_nullspace_C_calculation',r.FFT_nullspace_C_calculation, ...
    'sketched_SVD',r.sketched_SVD,'sketch_dim',r.sketch_dim, ...
    'visualize_C_matrix_sv',r.visualize_C_matrix_sv);

G=utils.G_matrices_3D(kCal,gridSize(1),gridSize(2),gridSize(3),r.tau,U, ...
    r.FFT_nullspace_C_calculation,'kernel_shape',r.kernel_shape, ...
    'FFT_interpolation',r.FFT_interpolation,'interp_zp',r.interp_zp, ...
    'sketched_SVD',r.sketched_SVD);
patchSize=size(U,1)/size(kCal,4);
clear U

[ST_maps,eigenValues,diagnostics]=nullspace_vectors_G_matrix_3D_STM_v2( ...
    kCal,dim_stm(1),dim_stm(2),dim_stm(3),L,G,patchSize, ...
    'Solver',r.Solver,'M',r.M,'Tolerance',r.Tolerance, ...
    'CheckEvery',r.CheckEvery,'RandomSeed',r.RandomSeed, ...
    'RayleighRitz',r.RayleighRitz,'FFT_interpolation',r.FFT_interpolation, ...
    'gauss_win_param',r.gauss_win_param,'PhaseAlignment',r.PhaseAlignment, ...
    'InterpolateEigenvalues',r.InterpolateEigenvalues, ...
    'NormalizeMaps',r.NormalizeMaps,'verbose',r.verbose);
clear G

if r.verbose, fprintf('Total 3D STM prototype time: %.3f s.\n',toc(totalTimer)); end
end