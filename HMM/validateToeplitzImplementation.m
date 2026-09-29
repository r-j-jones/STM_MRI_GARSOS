%% Part 1

radialOp = result.stmOperator.radialOperator;
t = 1;
N = radialOp.imageSize;
omega = double(radialOp.plans(t).omega);

% Replace these with the exact parameters from your plan builder.
nufftArguments = { ...
    N, [6 6], 2*N, N/2, ...
    'table', 2^10, 'minmax:kb'};

G = Gnufft(true(N), {omega, nufftArguments{:}});

x = complex( ...
    randn(N, 'single'), ...
    randn(N, 'single'));

% Existing low-level operation.
yLowLevel = nufft(x, radialOp.plans(t).st);

% Gnufft operation; full mask means x(:) is accepted.
yGnufft = G * x(:);

forwardError = norm(double(yLowLevel(:)) - double(yGnufft(:))) / ...
    max(norm(double(yLowLevel(:))), eps);

fprintf('Low-level NUFFT versus Gnufft error: %.3e\n', ...
    forwardError);



%% Part 2

normalLowLevel = nufft_adj( ...
    nufft(x, radialOp.plans(t).st), ...
    radialOp.plans(t).st);

normalGnufft = reshape(G' * (G * x(:)), N);

normalError = norm( ...
    double(normalLowLevel(:)) - double(normalGnufft(:))) / ...
    max(norm(double(normalLowLevel(:))), eps);

fprintf('Low-level versus Gnufft normal error: %.3e\n', ...
    normalError);


%% PART 3

tic;
T = build_gram(G, 1);
fprintf('One Gram build time: %.3f s\n', toc);

exact = reshape(G' * (G * x(:)), N);
toeplitzResult = reshape(T * x(:), N);

toeplitzError = norm( ...
    double(exact(:)) - double(toeplitzResult(:))) / ...
    max(norm(double(exact(:))), eps);

fprintf('Gnufft normal versus Toeplitz error: %.3e\n', ...
    toeplitzError);




%% PART 4

u = complex(randn(N), randn(N));
v = complex(randn(N), randn(N));

Tu = reshape(T * u(:), N);
Tv = reshape(T * v(:), N);

left = sum(conj(u(:)) .* Tv(:));
right = conj(sum(conj(v(:)) .* Tu(:)));

hermitianError = abs(left-right) / ...
    max([abs(left), abs(right), eps]);

energy = real(sum(conj(u(:)) .* Tu(:)));

fprintf('Toeplitz Hermitian error: %.3e\n', hermitianError);
fprintf('Toeplitz quadratic energy: %.6g\n', energy);



%% PART 5: 

radialOp = result.stmOperator.radialOperator;
stmOp = result.stmOperator;

toeplitzGram = buildGARSOSToeplitzGramOperators( ...
    radialOp, ...
    'Jd', [6 6], ...          % Replace with your actual values
    'Kd', 2 * radialOp.imageSize, ...
    'NShift', radialOp.imageSize / 2, ...
    'Kernel', 'minmax:kb', ...
    'Table', 2^10, ...
    'Verbose', true);

radialToeplitz = buildGARSOSRadialToeplitzNormalOperator( ...
    radialOp, toeplitzGram);

stmToeplitz = buildGARSOSSTMToeplitzNormalOperator( ...
    stmOp, radialToeplitz);

%% PART 6

nCoefficient = prod([stmOp.imageSize, stmOp.nBasis]);

c = complex( ...
    randn(nCoefficient, 1, 'single'), ...
    randn(nCoefficient, 1, 'single'));

tic;
exact = stmOp.normalVector(c);
exactTime = toc;

tic;
approximate = stmToeplitz.applyVector(c);
toeplitzTime = toc;

relativeError = norm(double(exact-approximate)) / ...
    max(norm(double(exact)), eps);

fprintf('Full STM Toeplitz relative error: %.3e\n', relativeError);
fprintf('Exact normal time: %.3f s\n', exactTime);
fprintf('Toeplitz normal time: %.3f s\n', toeplitzTime);
fprintf('Speedup: %.2fx\n', exactTime / toeplitzTime);



