function [ST_maps,eigenVal,diagnostics] = nullspace_vectors_G_matrix_3D_STM_v2( ...
    kCal,N1,N2,N3,L,G,patchSize,varargin)
%NULLSPACE_VECTORS_G_MATRIX_3D_STM_V2 Compute L 3D STMs from local G matrices.
%
% Implements:
%   Phase A: paged SVD reference.
%   Phase B: block orthogonal/subspace iteration.
%   Phase C: optional Rayleigh-Ritz refinement.
%   Phase D: optional PISCO-style 3D FFT interpolation and phase alignment.
%
% Inputs
%   kCal       N1_cal x N2_cal x N3_cal x Nc calibration data. In the STM
%              interpretation, dimension 4 contains virtual channels/time.
%   N1,N2,N3   Desired output spatial dimensions.
%   L          Number of STM vectors, 1 <= L <= Nc.
%   G          N1_g x N2_g x N3_g x Nc x Nc local matrices.
%   patchSize  Positive normalization scalar.
%
% Options
%   Solver                 'svd' (default) or 'orthogonal-iteration'.
%   M                      Maximum iteration count (default 100).
%   Tolerance              Maximum relative subspace residual (default 1e-8).
%   CheckEvery             Residual-check interval (default 5).
%   RandomSeed             Iterative initialization seed (default 0).
%   RayleighRitz           Refine/order the iterative basis (default true).
%   FFT_interpolation      Interpolate dimensions 1-3 (default true).
%   gauss_win_param        Gaussian-window parameter (default 100).
%   PhaseAlignment         'calibration' (default), 'first-channel', or 'none'.
%   InterpolateEigenvalues Interpolate selected eigenvalues (default true).
%   NormalizeMaps          Normalize each STM across Nc (default false).
%   verbose                Display iterative progress (default true).
%
% Outputs
%   ST_maps     N1 x N2 x N3 x Nc x L with interpolation, otherwise
%               N1_g x N2_g x N3_g x Nc x L.
%   eigenVal    Corresponding spatial grid x L, ordered ascending when SVD
%               or RayleighRitz is enabled.
%   diagnostics Solver information and low-grid residuals.

p = inputParser;
p.FunctionName = mfilename;
p.addRequired('kCal',@(x) isnumeric(x) && ndims(x)==4);
p.addRequired('N1',@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addRequired('N2',@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addRequired('N3',@(x) isnumeric(x) && isscalar(x) && x>=1);
p.addRequired('L',@(x) isnumeric(x) && isscalar(x) && x>=1 && x==floor(x));
p.addRequired('G',@(x) isnumeric(x) && ndims(x)==5);
p.addRequired('patchSize',@(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('Solver','svd',@(x) ischar(x)||(isstring(x)&&isscalar(x)));
p.addParameter('M',100,@(x) isnumeric(x)&&isscalar(x)&&x>=1&&x==floor(x));
p.addParameter('Tolerance',1e-8,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('CheckEvery',5,@(x) isnumeric(x)&&isscalar(x)&&x>=1&&x==floor(x));
p.addParameter('RandomSeed',0,@(x) isnumeric(x)&&isscalar(x));
p.addParameter('RayleighRitz',true,@(x) islogical(x)||(isnumeric(x)&&isscalar(x)));
p.addParameter('FFT_interpolation',true,@(x) islogical(x)||(isnumeric(x)&&isscalar(x)));
p.addParameter('gauss_win_param',100,@(x) isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('PhaseAlignment','calibration',@(x) ischar(x)||(isstring(x)&&isscalar(x)));
p.addParameter('InterpolateEigenvalues',true,@(x) islogical(x)||(isnumeric(x)&&isscalar(x)));
p.addParameter('NormalizeMaps',false,@(x) islogical(x)||(isnumeric(x)&&isscalar(x)));
p.addParameter('verbose',true,@(x) islogical(x)||(isnumeric(x)&&isscalar(x)));
p.parse(kCal,N1,N2,N3,L,G,patchSize,varargin{:});

solver = lower(char(p.Results.Solver));
if ~ismember(solver,{'svd','orthogonal-iteration'})
    error('Solver must be ''svd'' or ''orthogonal-iteration''.');
end
phaseMode = lower(char(p.Results.PhaseAlignment));
if ~ismember(phaseMode,{'calibration','first-channel','none'})
    error('PhaseAlignment must be ''calibration'', ''first-channel'', or ''none''.');
end

N1_g=size(G,1); N2_g=size(G,2); N3_g=size(G,3); Nc=size(G,4);
if size(G,5)~=Nc, error('G must end in equal Nc x Nc dimensions.'); end
if L>Nc, error('L must not exceed Nc.'); end
if size(kCal,4)~=Nc, error('size(kCal,4) must equal Nc.'); end
if ~logical(p.Results.FFT_interpolation) && any([N1 N2 N3]~=[N1_g N2_g N3_g])
    error('Without FFT interpolation, requested dimensions must equal the G grid.');
end

numVoxels=N1_g*N2_g*N3_g;
A=reshape(permute(G,[4 5 1 2 3]),[Nc Nc numVoxels])/patchSize;
A=0.5*(A+pagectranspose_local(A));
clear G

diagnostics=struct('solver',solver,'iterations',0,'converged',true, ...
    'rayleighRitz',logical(p.Results.RayleighRitz), ...
    'relativeResidual',zeros(numVoxels,1));

switch solver
    case 'svd'
        [~,s,V]=pagesvd(A,'econ','vector');
        idx=Nc:-1:Nc-L+1;
        S=V(:,idx,:);
        eigenPages=reshape(s(idx,:,:),[L numVoxels]);
        diagnostics.relativeResidual=subspace_residual(A,S);

    case 'orthogonal-iteration'
        prior=rng; restoreRng=onCleanup(@()rng(prior)); %#ok<NASGU>
        rng(p.Results.RandomSeed,'twister');
        S=complex(randn(Nc,L,numVoxels,'like',A),randn(Nc,L,numVoxels,'like',A));
        S=orthonormalize_pages(S);
        alpha=page_max_eigenvalue(A);
        B=repmat(eye(Nc,'like',A),[1 1 numVoxels]).*reshape(alpha,[1 1 numVoxels])-A;
        converged=false;
        for iter=1:p.Results.M
            S=orthonormalize_pages(pagemtimes(B,S));
            if p.Results.Tolerance>0 && (mod(iter,p.Results.CheckEvery)==0 || iter==p.Results.M)
                residual=subspace_residual(A,S);
                maxResidual=max(gather_if_needed(residual));
                if p.Results.verbose
                    fprintf('3D STM orthogonal iteration %d: max residual %.3e\n',iter,maxResidual);
                end
                if maxResidual<=p.Results.Tolerance, converged=true; break; end
            end
        end
        H=pagemtimes(S,'ctranspose',pagemtimes(A,S),'none');
        H=0.5*(H+pagectranspose_local(H));
        if logical(p.Results.RayleighRitz)
            [Uritz,eigenPages]=pageeig_small(H);
            S=pagemtimes(S,Uritz);
        else
            eigenPages=reshape(real(pagediagonal(H)),[L numVoxels]);
        end
        diagnostics.iterations=iter;
        diagnostics.converged=converged || p.Results.Tolerance==0;
        diagnostics.relativeResidual=subspace_residual(A,S);
end

ST_maps=reshape(permute(S,[3 1 2]),[N1_g N2_g N3_g Nc L]);
eigenVal=reshape(permute(eigenPages,[2 1]),[N1_g N2_g N3_g L]);

if logical(p.Results.FFT_interpolation)
    [ST_maps,eigenVal]=interpolate_3d_stm(ST_maps,eigenVal,kCal, ...
        [N1 N2 N3],p.Results.gauss_win_param,phaseMode, ...
        logical(p.Results.InterpolateEigenvalues));
else
    ST_maps=align_phase(ST_maps,kCal,phaseMode);
end

if logical(p.Results.NormalizeMaps)
    den=sqrt(sum(abs(ST_maps).^2,4));
    den(den==0)=1;
    ST_maps=ST_maps./den;
end
end

function [mapsOut,eigOut]=interpolate_3d_stm(maps,eigVal,kCal,outSize,gaussParam,phaseMode,interpEig)
[N1g,N2g,N3g,Nc,L]=size(maps);
[N1c,N2c,N3c,~]=size(kCal);
wx=hamming_local(N1g); wy=hamming_local(N2g); wz=hamming_local(N3g);
w=reshape(wx,[N1g 1 1]).*reshape(wy,[1 N2g 1]).*reshape(wz,[1 1 N3g]);

% Align low-grid eigenvectors before Fourier interpolation.
maps=align_phase(maps,kCal,phaseMode);
winMaps=reshape(w,[N1g N2g N3g 1 1]);
mapsOut=ifft3_centered_zp(fft3_centered(maps).*winMaps,outSize);

if interpEig
    eigOut=abs(ifft3_centered_zp(fft3_centered(eigVal).*reshape(w,[N1g N2g N3g 1]),outSize));
    scale=max(eigOut,[],[1 2 3]); scale(scale==0)=1;
    eigOut=eigOut./scale;
else
    eigOut=eigVal;
end

% Nested helper data are validated here to catch unsupported grids early.
if any([N1c N2c N3c]>[N1g N2g N3g])
    error('Calibration dimensions must not exceed the low-resolution G grid.');
end
end

function maps=align_phase(maps,kCal,mode)
if strcmp(mode,'none'), return; end
[N1g,N2g,N3g,Nc,L]=size(maps);
if strcmp(mode,'first-channel')
    ref=exp(-1i*angle(maps(:,:,:,1,:)));
    maps=maps.*ref;
    return
end
[N1c,N2c,N3c,~]=size(kCal);
if any([N1c N2c N3c]>[N1g N2g N3g])
    error('Calibration dimensions must not exceed the map grid for calibration phase alignment.');
end
im=zeros(N1g,N2g,N3g,Nc,'like',maps);
ix=center_indices(N1g,N1c); iy=center_indices(N2g,N2c); iz=center_indices(N3g,N3c);
im(ix,iy,iz,:)=kCal;
win=reshape(gausswin_local(N1g,100),[N1g 1 1]).* ...
    reshape(gausswin_local(N2g,100),[1 N2g 1]).* ...
    reshape(gausswin_local(N3g,100),[1 1 N3g]);
im=ifft3_centered(im.*win);
im=reshape(im,[N1g N2g N3g Nc 1]);
num=sum(conj(maps).*im,4); den=sum(abs(maps).^2,4); den(den==0)=1;
maps=maps.*exp(1i*angle(num./den));
end

function idx=center_indices(big,small)
start=floor((big-small)/2)+1; idx=start:start+small-1;
end
function y=fft3_centered(x), y=fftshift(fftn(ifftshift(x,1),[],1),1); y=fftshift(fftn(ifftshift(y,2),[],2),2); y=fftshift(fftn(ifftshift(y,3),[],3),3); end
function y=ifft3_centered(x), y=fftshift(ifftn(ifftshift(x,1),[],1),1); y=fftshift(ifftn(ifftshift(y,2),[],2),2); y=fftshift(ifftn(ifftshift(y,3),[],3),3); end
function y=ifft3_centered_zp(x,outSize)
y=fftshift(ifftn(ifftshift(x,1),outSize(1),1),1); y=fftshift(ifftn(ifftshift(y,2),outSize(2),2),2); y=fftshift(ifftn(ifftshift(y,3),outSize(3),3),3);
end
function w=hamming_local(n), if n==1, w=1; else, q=(0:n-1)'; w=0.54-0.46*cos(2*pi*q/(n-1)); end, end
function w=gausswin_local(n,a), if n==1, w=1; else, q=linspace(-1,1,n)'; w=exp(-0.5*(a*q).^2); end, end
function Q=orthonormalize_pages(X)
[~,L,~]=size(X); Q=X;
for pass=1:2
 for l=1:L
  for k=1:l-1
   c=pagemtimes(Q(:,k,:),'ctranspose',Q(:,l,:),'none'); Q(:,l,:)=Q(:,l,:)-Q(:,k,:).*c;
  end
  n=sqrt(real(pagemtimes(Q(:,l,:),'ctranspose',Q(:,l,:),'none'))); n=max(n,eps(class(Q))); Q(:,l,:)=Q(:,l,:)./n;
 end
end
end
function r=subspace_residual(A,S)
AS=pagemtimes(A,S); H=pagemtimes(S,'ctranspose',AS,'none'); R=AS-pagemtimes(S,H);
r=reshape(sqrt(sum(abs(R).^2,[1 2]))./max(sqrt(sum(abs(AS).^2,[1 2])),eps(class(A))),[],1);
end
function alpha=page_max_eigenvalue(A)
np=size(A,3); alpha=zeros(np,1,'like',real(A));
for q=1:np, alpha(q)=max(real(eig(A(:,:,q),'vector'))); end
alpha=alpha+eps(class(A)).*max(abs(alpha),1);
end
function [U,D]=pageeig_small(H)
L=size(H,1); np=size(H,3); U=zeros(L,L,np,'like',H); D=zeros(L,np,'like',real(H));
for q=1:np, [Uq,dq]=eig(H(:,:,q),'vector'); [dq,o]=sort(real(dq),'ascend'); U(:,:,q)=Uq(:,o); D(:,q)=dq; end
end
function d=pagediagonal(H)
L=size(H,1); np=size(H,3); d=zeros(L,1,np,'like',real(H)); for l=1:L, d(l,1,:)=real(H(l,l,:)); end
end
function XH=pagectranspose_local(X), XH=permute(conj(X),[2 1 3]); end
function X=gather_if_needed(X), if isa(X,'gpuArray'), X=gather(X); end, end