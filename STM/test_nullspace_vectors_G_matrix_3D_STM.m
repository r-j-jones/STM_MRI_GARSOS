function test_nullspace_vectors_G_matrix_3D_STM
%TEST_NULLSPACE_VECTORS_G_MATRIX_3D_STM Synthetic Phase A/B validation.
%
% Place this test and nullspace_vectors_G_matrix_3D_STM.m on the MATLAB
% path, then run:
%   test_nullspace_vectors_G_matrix_3D_STM

rng(7,'twister');
N1_g = 3;
N2_g = 2;
N3_g = 2;
Nc = 5;
L = 2;
patchSize = 1;
numVoxels = N1_g*N2_g*N3_g;

% Build pagewise Hermitian positive-semidefinite matrices with a known,
% separated spectrum. The two smallest eigenvalues define the target STM
% subspace at each voxel.
Gpages = zeros(Nc,Nc,numVoxels);
trueValues = zeros(Nc,numVoxels);
for q = 1:numVoxels
    [Q,~] = qr(randn(Nc) + 1i*randn(Nc),0);
    d = [0.01 + 0.001*q; 0.04 + 0.001*q; 0.35; 0.70; 1.10];
    Gpages(:,:,q) = Q*diag(d)*Q';
    trueValues(:,q) = d;
end
G = reshape(permute(Gpages,[3 1 2]),[N1_g N2_g N3_g Nc Nc]);

% kCal is not used when FFT interpolation is disabled, but its fourth
% dimension must match Nc to enforce the intended interface.
kCal = complex(zeros(2,2,2,Nc));

[ST_svd,eig_svd,diag_svd] = nullspace_vectors_G_matrix_3D_STM( ...
    kCal,N1_g,N2_g,N3_g,L,G,patchSize, ...
    'Solver','svd','FFT_interpolation',0,'verbose',0);

[ST_oi,eig_oi,diag_oi] = nullspace_vectors_G_matrix_3D_STM( ...
    kCal,N1_g,N2_g,N3_g,L,G,patchSize, ...
    'Solver','orthogonal-iteration','M',300,'Tolerance',1e-11, ...
    'CheckEvery',5,'RandomSeed',11,'FFT_interpolation',0,'verbose',0);

assert(isequal(size(ST_svd),[N1_g N2_g N3_g Nc L]), ...
    'Unexpected Phase A map dimensions.');
assert(isequal(size(ST_oi),[N1_g N2_g N3_g Nc L]), ...
    'Unexpected Phase B map dimensions.');
assert(isequal(size(eig_svd),[N1_g N2_g N3_g L]), ...
    'Unexpected Phase A eigenvalue dimensions.');
assert(isequal(size(eig_oi),[N1_g N2_g N3_g L]), ...
    'Unexpected Phase B eigenvalue dimensions.');
assert(all(diff(reshape(eig_svd,[],L),1,2) >= -1e-12,'all'), ...
    'Phase A eigenvalues are not ascending.');
assert(all(diff(reshape(eig_oi,[],L),1,2) >= -1e-10,'all'), ...
    'Phase B eigenvalues are not ascending.');

Ssvd = maps_to_pages(ST_svd);
Soi = maps_to_pages(ST_oi);
maxOrthSVD = 0;
maxOrthOI = 0;
maxResidualSVD = 0;
maxResidualOI = 0;
maxProjectorDifference = 0;
maxEigAgreement = 0;

for q = 1:numVoxels
    A = Gpages(:,:,q)/patchSize;
    Vs = Ssvd(:,:,q);
    Vo = Soi(:,:,q);
    ls = reshape(eig_svd(q),1,[]); %#ok<NASGU>

    maxOrthSVD = max(maxOrthSVD,norm(Vs'*Vs-eye(L),'fro'));
    maxOrthOI = max(maxOrthOI,norm(Vo'*Vo-eye(L),'fro'));

    lambdaS = diag(Vs'*A*Vs);
    lambdaO = diag(Vo'*A*Vo);
    for l = 1:L
        maxResidualSVD = max(maxResidualSVD, ...
            norm(A*Vs(:,l)-lambdaS(l)*Vs(:,l))/max(norm(A*Vs(:,l)),eps));
        maxResidualOI = max(maxResidualOI, ...
            norm(A*Vo(:,l)-lambdaO(l)*Vo(:,l))/max(norm(A*Vo(:,l)),eps));
    end

    maxProjectorDifference = max(maxProjectorDifference, ...
        norm(Vs*Vs'-Vo*Vo','fro'));

    exact = sort(real(eig(A,'vector')),'ascend');
    maxEigAgreement = max(maxEigAgreement, ...
        max(abs(sort(real(lambdaO),'ascend')-exact(1:L))));
end

assert(maxOrthSVD < 1e-10,'Phase A vectors are not orthonormal.');
assert(maxOrthOI < 1e-9,'Phase B vectors are not orthonormal.');
assert(maxResidualSVD < 1e-9,'Phase A eigen residual is too large.');
assert(maxResidualOI < 1e-7,'Phase B eigen residual is too large.');
assert(maxProjectorDifference < 1e-7, ...
    'Phase A/B invariant subspaces do not agree.');
assert(maxEigAgreement < 1e-8, ...
    'Phase B Ritz values do not agree with MATLAB eig.');
assert(max(abs(reshape(eig_svd,[],L)-trueValues(1:L,:).')) < 1e-10, ...
    'Phase A selected eigenvalues do not match the known spectrum.');

% L=1 compatibility and dimensions.
[ST_one,eig_one] = nullspace_vectors_G_matrix_3D_STM( ...
    kCal,N1_g,N2_g,N3_g,1,G,patchSize, ...
    'Solver','svd','FFT_interpolation',0,'verbose',0);
assert(isequal(size(ST_one),[N1_g N2_g N3_g Nc 1]), ...
    'L=1 map dimensions are incorrect.');
assert(isequal(size(eig_one),[N1_g N2_g N3_g 1]), ...
    'L=1 eigenvalue dimensions are incorrect.');

fprintf('3D STM Phase A/B synthetic test passed.\n');
fprintf('  SVD max orthonormality error: %.3e\n',maxOrthSVD);
fprintf('  OI  max orthonormality error: %.3e\n',maxOrthOI);
fprintf('  SVD max eigen residual:       %.3e\n',maxResidualSVD);
fprintf('  OI  max eigen residual:       %.3e\n',maxResidualOI);
fprintf('  Max projector difference:     %.3e\n',maxProjectorDifference);
fprintf('  OI iterations: %d; converged: %d\n', ...
    diag_oi.iterations,diag_oi.converged);
fprintf('  SVD max subspace residual:    %.3e\n', ...
    max(diag_svd.relativeResidual));
end

function S = maps_to_pages(ST_maps)
sz = size(ST_maps);
N1 = sz(1);
N2 = sz(2);
N3 = sz(3);
Nc = sz(4);
L = sz(5);
S = permute(reshape(ST_maps,[N1*N2*N3 Nc L]),[2 3 1]);
end