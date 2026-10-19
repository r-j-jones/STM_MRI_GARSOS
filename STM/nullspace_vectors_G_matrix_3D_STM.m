function [ST_maps,eigenVal,diagnostics] = nullspace_vectors_G_matrix_3D_STM( ...
    kCal,N1,N2,N3,L,G,patchSize,varargin)
%NULLSPACE_VECTORS_G_MATRIX_3D_STM Compute L low-eigenvalue 3D STM vectors.
%
% Phase A provides a paged-SVD reference. Phase B provides block
% orthogonal/subspace iteration followed by a Rayleigh-Ritz refinement.
% FFT interpolation is intentionally deferred; this function currently
% returns maps on the spatial grid represented by G.
%
% Inputs
%   kCal       N1_cal x N2_cal x N3_cal x Nc calibration array. For STM,
%              dimension 4 represents virtual channels/time frames.
%   N1,N2,N3   Requested full-resolution dimensions. Retained for interface
%              compatibility; Phase A/B require them to equal the G grid
%              whenever FFT_interpolation is zero.
%   L          Number of STMs, integer in [1,Nc].
%   G          N1_g x N2_g x N3_g x Nc x Nc local G matrices.
%   patchSize  Positive normalization scalar used by PISCO/STM.
%
% Name-value options
%   'Solver'             'svd' (default) or 'orthogonal-iteration'.
%   'M'                  Maximum block-iteration count (default 100).
%   'Tolerance'          Relative maximum subspace residual tolerance
%                        (default 1e-8). Set to 0 to force M iterations.
%   'CheckEvery'         Residual-check interval (default 5).
%   'RandomSeed'         Seed used for iterative initialization (default 0).
%   'FFT_interpolation'  Must be 0 for Phase A/B (default 0).
%   'verbose'            Display iterative progress (default 1).
%
% Outputs
%   ST_maps     N1_g x N2_g x N3_g x Nc x L, with L ordered by ascending
%               eigenvalue of G/patchSize.
%   eigenVal    N1_g x N2_g x N3_g x L selected normalized eigenvalues.
%   diagnostics Structure containing solver, iteration count, convergence,
%               and final relative subspace residuals.
%
% Notes
%   Individual complex eigenvectors are defined only up to phase. For
%   repeated/clustered eigenvalues, compare invariant-subspace projectors
%   rather than individual vectors.

p = inputParser;
p.FunctionName = mfilename;
p.addRequired('kCal',@(x) isnumeric(x) && ndims(x) == 4);
p.addRequired('N1',@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
p.addRequired('N2',@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
p.addRequired('N3',@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
p.addRequired('L',@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
p.addRequired('G',@(x) isnumeric(x) && ndims(x) == 5);
p.addRequired('patchSize',@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
p.addParameter('Solver','svd',@(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('M',100,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
p.addParameter('Tolerance',1e-8,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('CheckEvery',5,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
p.addParameter('RandomSeed',0,@(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('FFT_interpolation',0,@(x) isnumeric(x) || islogical(x));
p.addParameter('verbose',1,@(x) isnumeric(x) || islogical(x));
p.parse(kCal,N1,N2,N3,L,G,patchSize,varargin{:});

solver = lower(char(p.Results.Solver));
validSolvers = {'svd','orthogonal-iteration'};
if ~ismember(solver,validSolvers)
    error('Unsupported Solver ''%s''. Use ''svd'' or ''orthogonal-iteration''.',solver);
end
if logical(p.Results.FFT_interpolation)
    error(['FFT interpolation is not implemented in Phase A/B. Set ' ...
        '''FFT_interpolation'',0 and validate the low-resolution eigensolver first.']);
end

N1_g = size(G,1);
N2_g = size(G,2);
N3_g = size(G,3);
Nc = size(G,4);
if size(G,5) ~= Nc
    error('G must have equal channel dimensions: [...,Nc,Nc].');
end
if L > Nc
    error('L (%d) must not exceed Nc (%d).',L,Nc);
end
if size(kCal,4) ~= Nc
    error('size(kCal,4) (%d) must equal the G channel dimension Nc (%d).',size(kCal,4),Nc);
end
if any([N1 N2 N3] ~= [N1_g N2_g N3_g])
    error(['Without FFT interpolation, [N1 N2 N3] must equal the G grid ' ...
        '[N1_g N2_g N3_g].']);
end

numVoxels = N1_g*N2_g*N3_g;
Gpages = reshape(permute(G,[4 5 1 2 3]),[Nc Nc numVoxels]);
% Suppress small numerical departures from Hermitian symmetry.
Gpages = 0.5*(Gpages + pagectranspose_local(Gpages));
A = Gpages/p.Results.patchSize;
clear Gpages

diagnostics = struct('solver',solver,'iterations',0,'converged',true, ...
    'relativeResidual',zeros(numVoxels,1));

switch solver
    case 'svd'
        [~,singularValues,Vpage] = pagesvd(A,'econ','vector');
        idx = Nc:-1:(Nc-L+1); % ascending selected eigenvalue order
        S = Vpage(:,idx,:);
        selectedValues = singularValues(idx,:,:);
        eigenPages = reshape(selectedValues,[L numVoxels]);
        diagnostics.relativeResidual = subspace_residual(A,S);

    case 'orthogonal-iteration'
        oldRng = rng;
        cleanup = onCleanup(@() rng(oldRng)); %#ok<NASGU>
        rng(p.Results.RandomSeed,'twister');

        if isa(A,'gpuArray')
            baseClass = classUnderlying(A);
            S = complex(gpuArray.randn(Nc,L,numVoxels,baseClass), ...
                gpuArray.randn(Nc,L,numVoxels,baseClass));
        else
            S = complex(randn(Nc,L,numVoxels,'like',A), ...
                randn(Nc,L,numVoxels,'like',A));
        end
        S = orthonormalize_pages(S);

        % A safe per-page shift. B and A have the same eigenvectors, while
        % the smallest eigenvalues of A become the largest eigenvalues of B.
        alpha = page_max_eigenvalue(A);
        Ipages = repmat(eye(Nc,'like',A),[1 1 numVoxels]);
        B = Ipages.*reshape(alpha,[1 1 numVoxels]) - A;
        clear Ipages

        converged = false;
        residual = inf(numVoxels,1);
        for iter = 1:p.Results.M
            S = pagemtimes(B,S);
            S = orthonormalize_pages(S);

            shouldCheck = p.Results.Tolerance > 0 && ...
                (mod(iter,p.Results.CheckEvery) == 0 || iter == p.Results.M);
            if shouldCheck
                residual = subspace_residual(A,S);
                maxResidual = max(gather_if_needed(residual));
                if logical(p.Results.verbose)
                    fprintf('3D STM orthogonal iteration %d: max relative residual %.3e\n', ...
                        iter,maxResidual);
                end
                if maxResidual <= p.Results.Tolerance
                    converged = true;
                    break
                end
            end
        end

        % Rayleigh-Ritz: rotate the converged subspace into ordered
        % approximate eigenvectors and obtain their Ritz values.
        H = pagemtimes(S,'ctranspose',pagemtimes(A,S),'none');
        H = 0.5*(H + pagectranspose_local(H));
        [Uritz,Dritz] = pageeig_small(H);
        S = pagemtimes(S,Uritz);
        eigenPages = Dritz;
        residual = subspace_residual(A,S);

        diagnostics.iterations = iter;
        diagnostics.converged = converged || p.Results.Tolerance == 0;
        diagnostics.relativeResidual = residual;
        if p.Results.Tolerance > 0 && ~converged && logical(p.Results.verbose)
            warning('nullspace_vectors_G_matrix_3D_STM:NotConverged', ...
                ['Maximum iteration count reached. Final maximum relative ' ...
                 'subspace residual is %.3e.'],max(gather_if_needed(residual)));
        end
end

ST_maps = reshape(permute(S,[3 1 2]),[N1_g N2_g N3_g Nc L]);
eigenVal = reshape(permute(eigenPages,[2 1]),[N1_g N2_g N3_g L]);
end

function Q = orthonormalize_pages(X)
% Modified Gram-Schmidt with a second pass for improved stability.
[~,L,~] = size(X);
Q = X;
for pass = 1:2
    for l = 1:L
        for k = 1:l-1
            coeff = pagemtimes(Q(:,k,:),'ctranspose',Q(:,l,:),'none');
            Q(:,l,:) = Q(:,l,:) - Q(:,k,:).*coeff;
        end
        nrm = sqrt(real(pagemtimes(Q(:,l,:),'ctranspose',Q(:,l,:),'none')));
        floorValue = eps(class_underlying_local(Q));
        nrm = max(nrm,floorValue);
        Q(:,l,:) = Q(:,l,:)./nrm;
    end
end
end

function residual = subspace_residual(A,S)
AS = pagemtimes(A,S);
H = pagemtimes(S,'ctranspose',AS,'none');
R = AS - pagemtimes(S,H);
num = sqrt(sum(abs(R).^2,[1 2]));
den = max(sqrt(sum(abs(AS).^2,[1 2])),eps(class_underlying_local(A)));
residual = reshape(num./den,[],1);
end

function alpha = page_max_eigenvalue(A)
% Exact pagewise maximum eigenvalue gives a robust shift for Phase B.
numPages = size(A,3);
alpha = zeros(numPages,1,'like',real(A));
for q = 1:numPages
    values = eig(A(:,:,q),'vector');
    alpha(q) = max(real(values));
end
% A small positive pad avoids a zero transformed leading eigenvalue.
pad = eps(class_underlying_local(A)).*max(abs(alpha),1);
alpha = alpha + pad;
end

function [U,D] = pageeig_small(H)
% Pagewise eigendecomposition of small L-by-L Hermitian matrices.
L = size(H,1);
numPages = size(H,3);
U = zeros(L,L,numPages,'like',H);
D = zeros(L,numPages,'like',real(H));
for q = 1:numPages
    [Uq,Dq] = eig(H(:,:,q),'vector');
    [dq,order] = sort(real(Dq),'ascend');
    U(:,:,q) = Uq(:,order);
    D(:,q) = dq;
end
end

function XH = pagectranspose_local(X)
XH = permute(conj(X),[2 1 3]);
end

function name = class_underlying_local(X)
if isa(X,'gpuArray')
    name = classUnderlying(X);
else
    name = class(X);
end
end

function X = gather_if_needed(X)
if isa(X,'gpuArray')
    X = gather(X);
end
end