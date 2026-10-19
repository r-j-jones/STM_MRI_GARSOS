function test_nullspace_vectors_G_matrix_3D_STM_v2
% Synthetic tests for Phases A-D.
rng(9);
lo=[3 3 2]; hi=[5 4 4]; Nc=5; L=2; np=prod(lo);
Gp=zeros(Nc,Nc,np); known=zeros(np,L);
for q=1:np
 [Q,~]=qr(randn(Nc)+1i*randn(Nc),0);
 d=[.01+.0001*q,.05+.0001*q,.3,.7,1.1];
 Gp(:,:,q)=Q*diag(d)*Q'; known(q,:)=d(1:L);
end
G=reshape(permute(Gp,[3 1 2]),[lo Nc Nc]);
kCal=complex(randn(2,2,2,Nc),randn(2,2,2,Nc));

[A,ea]=nullspace_vectors_G_matrix_3D_STM_v2(kCal,lo(1),lo(2),lo(3),L,G,1, ...
 'Solver','svd','FFT_interpolation',0,'PhaseAlignment','none');
[B,eb,db]=nullspace_vectors_G_matrix_3D_STM_v2(kCal,lo(1),lo(2),lo(3),L,G,1, ...
 'Solver','orthogonal-iteration','M',300,'Tolerance',1e-10, ...
 'RayleighRitz',1,'FFT_interpolation',0,'PhaseAlignment','none','verbose',0);
assert(isequal(size(A),[lo Nc L])); assert(isequal(size(B),[lo Nc L]));
assert(max(abs(reshape(ea,[],L)-known),[],'all')<1e-10);
assert(all(diff(reshape(eb,[],L),1,2)>=-1e-10,'all'));

Ap=reshape(permute(reshape(A,[np Nc L]),[2 3 1]),[Nc L np]);
Bp=reshape(permute(reshape(B,[np Nc L]),[2 3 1]),[Nc L np]);
maxProj=0;
for q=1:np, maxProj=max(maxProj,norm(Ap(:,:,q)*Ap(:,:,q)'-Bp(:,:,q)*Bp(:,:,q)','fro')); end
assert(maxProj<1e-6); assert(db.converged);

% Rayleigh-Ritz-off option retains the target subspace but does not promise
% ordered individual vectors.
[C,~,dc]=nullspace_vectors_G_matrix_3D_STM_v2(kCal,lo(1),lo(2),lo(3),L,G,1, ...
 'Solver','orthogonal-iteration','M',300,'Tolerance',1e-10, ...
 'RayleighRitz',0,'FFT_interpolation',0,'PhaseAlignment','none','verbose',0);
assert(isequal(size(C),[lo Nc L])); assert(~dc.rayleighRitz);

% Phase D dimensions with interpolation; use first-channel alignment to keep
% this synthetic test independent of calibration-grid embedding semantics.
[D,ed]=nullspace_vectors_G_matrix_3D_STM_v2(kCal,hi(1),hi(2),hi(3),L,G,1, ...
 'Solver','svd','FFT_interpolation',1,'PhaseAlignment','first-channel', ...
 'InterpolateEigenvalues',1,'NormalizeMaps',1);
assert(isequal(size(D),[hi Nc L])); assert(isequal(size(ed),[hi L]));
assert(all(isfinite(D),'all')); assert(all(isfinite(ed),'all'));
mapNorm=sqrt(sum(abs(D).^2,4));
assert(max(abs(mapNorm(mapNorm>0)-1))<1e-8);

fprintf('3D STM v2 Phase A-D tests passed; max projector difference %.3e.\n',maxProj);
end