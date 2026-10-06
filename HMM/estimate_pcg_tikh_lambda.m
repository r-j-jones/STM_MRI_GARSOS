function tikh_lam = estimate_pcg_tikh_lambda( result, lamscale, niters, verbose )
%tikh_lam = estimate_pcg_tikh_lambda( result, lamscale, niters )
%
% result contains field stmOperator.normalVector
% lamcsale is desired lambda (def. 1e-2)
% niters is # of power iters
% 
% tikh_lam = spectralScale*lamscale;

if nargin<2
    lamscale = 1e-2;
end
if nargin<3
    niters = 10;
end
if nargin<4, verbose = true; end

if verbose, fprintf('Estimating PCG Tikh lambda..\n'); end


optest = result.stmOperator;
% normalOperator = @(x) optest.normalVector(x);
AHA = optest.normalVector;

nUnknown = prod(optest.imageSize) * optest.nBasis;
x = randn(nUnknown, 1, optest.precisionClass) + ...
    1i * randn(nUnknown, 1, optest.precisionClass);
x = x / norm(x);

estimates = zeros(15, 1);
for k = 1:niters
    if verbose
        fprintf(' iter %d/%d\t',k,niters)
        tic
    end
    y = AHA(x);
    if verbose, toc, end
    if verbose, tic, end
    
    estimates(k) = real(x' * y) / real(x' * x);
    x = y / norm(y);
    if verbose, toc, end
    if verbose, fprintf('%g\n',(estimates(k))); end
end

spectralScale = estimates(end);
tikh_lam = spectralScale * lamscale;
if verbose
    fprintf('Estimated ||A''A||_2  = %.6g\n', spectralScale);
    fprintf('Estimated Tikh lambda = %.6g\n', tikh_lam);
end
if verbose
    figure;
    plot(1:niters, estimates, '-o');
    ylabel('Est. spectral scale');
    xlabel('Iter');
    drawnow;
end

end


