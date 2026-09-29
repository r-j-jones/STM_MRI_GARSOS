function reconstruction = STM_reconstruction_GARSOS(result, opts)
%STM_RECONSTRUCTION_GARSOS Reconstruct GAR-SOS radial data with STM maps.
%
%   reconstruction = STM_reconstruction_GARSOS(result)
%   reconstruction = STM_reconstruction_GARSOS(result, opts)
%
% RESULT is the structure returned by main_STM_GARSOS.  In particular, it
% must contain the radial frame data and the composed STM/radial operator
% returned by runGARSOSSTMIntegration.  The reconstruction solves the
% multichannel radial least-squares problem in the STM coefficient domain,
% followed by an optional Tikhonov-regularized solve.
%
% Optional fields in OPTS:
%   tolerance       PCG stopping tolerance (default: 1e-6)
%   maxIterations   Maximum PCG iterations (default: 50)
%   lambda          Tikhonov parameter (default: 1e-3)
%   runTikhonov    Run the regularized reconstruction (default: true)
%   verbose         Display reconstruction progress messages (default: true)
%   initialGuess    Initial STM coefficient vector (default: zeros)
%   usePreconditioner Use STM-domain diagonal preconditioning. Default: false.
%   preconditionerNUFFTApproximation   'scalar' or 'image'. Default: 'scalar'.
%   preconditionerFloorFraction    Relative floor for the diagonal. Default: 1e-4.

    narginchk(1, 2);
    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    if ~isstruct(result)
        error('STM_reconstruction_GARSOS:InvalidResult', ...
            'result must be the structure returned by main_STM_GARSOS.');
    end
    if ~isstruct(opts)
        error('STM_reconstruction_GARSOS:InvalidOptions', ...
            'opts must be a structure.');
    end

    requiredFields = {'stmOperator', 'radialData', 'ST_maps'};
    if ~all(isfield(result, requiredFields))
        error('STM_reconstruction_GARSOS:MissingResultField', ...
            'result must contain stmOperator, radialData, and ST_maps.');
    end

    tolerance = getOption(opts, 'tolerance', 1e-5);
    maxIterations = getOption(opts, 'maxIterations', 50);
    lambda = getOption(opts, 'lambda', 1e-3);
    runTikhonov = getOption(opts, 'runTikhonov', true);
    verbose = getOption(opts, 'verbose', true);
    usePreconditioner = getOption(opts, ...
    'usePreconditioner', false);
    preconditionerApproximation = getOption(opts, ...
        'preconditionerNUFFTApproximation', 'scalar');
    preconditionerFloorFraction = getOption(opts, ...
        'preconditionerFloorFraction', 1e-4);

    validateattributes(tolerance, {'numeric'}, {'scalar', 'real', 'positive'});
    validateattributes(maxIterations, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(lambda, {'numeric'}, {'scalar', 'real', 'nonnegative'});
    validateattributes(runTikhonov, {'logical', 'numeric'}, {'scalar'});
    validateattributes(verbose, {'logical', 'numeric'}, {'scalar'});
    validateattributes(usePreconditioner, ...
        {'logical', 'numeric'}, {'scalar'});    
    validateattributes(preconditionerFloorFraction, {'numeric'}, ...
        {'scalar', 'real', 'positive'});


    op = result.stmOperator;
    data = result.radialData;  
    if isfield(op, 'useSinglePrecision') && op.useSinglePrecision
        data = single(data);
    end
    if ~isfield(op, 'forwardVector') || ~isfield(op, 'adjointVector') || ...
            ~isfield(op, 'normalVector')
        error('STM_reconstruction_GARSOS:InvalidOperator', ...
            'stmOperator must expose forwardVector, adjointVector, and normalVector.');
    end
    if ~isnumeric(data) || ~isequal(size(data), ...
            [op.radialOperator.nReadout, op.radialOperator.nSpokes, ...
             op.radialOperator.nCoil, op.nFrame])
        error('STM_reconstruction_GARSOS:DataSizeMismatch', ...
            'radialData dimensions do not match the radial encoding operator.');
    end

    coefficientSize = [size(result.ST_maps, 1), size(result.ST_maps, 2), ...
        size(result.ST_maps, 4)];
    nCoefficient = prod(coefficientSize);
    if isfield(opts, 'initialGuess') && ~isempty(opts.initialGuess)
        initialGuess = opts.initialGuess(:);
        if numel(initialGuess) ~= nCoefficient
            error('STM_reconstruction_GARSOS:InitialGuessSize', ...
                'initialGuess must contain one value per STM coefficient.');
        end
    else
        initialGuess = zeros(nCoefficient, 1, 'like', data);
    end

    rhs = op.adjointVector(data(:));
    if isfield(op, 'useSinglePrecision') && op.useSinglePrecision
        rhs = single(rhs);
    end

    if useToeplitzNormal
        normalOperator = stmToeplitz.applyVector;
    else
        normalOperator = @(x) op.normalVector(x);
    end

    if logical(verbose)
        disp('Starting CG-STM GAR-SOS reconstruction...');
    end
    tRecon = tic;
    [coefficients, flag, relativeResidual, iterations, residualHistory] = ...
        pcg(normalOperator, rhs, tolerance, maxIterations, [], [], initialGuess);
    if ~isscalar(flag)
        error('STM_reconstruction_GARSOS:PCGFailure', ...
            'The unregularized PCG solve returned an invalid status.');
    end

    reconstruction = makeOutput(coefficients, flag, relativeResidual, ...
        iterations, residualHistory, toc(tRecon), 'cg');

    if logical(runTikhonov)
        if logical(verbose)
            disp('Starting CG-STM GAR-SOS + Tikhonov reconstruction...');
        end
        tTik = tic;
        lambdaTyped = cast(lambda, 'like', rhs);
        tikOperator = @(x) normalOperator(x) + lambdaTyped .* x;
        % tikOperator = @(x) normalOperator(x) + lambda * x;
        [coefficientsTik, flagTik, relativeResidualTik, iterationsTik, ...
            residualHistoryTik] = pcg(tikOperator, rhs, tolerance, ...
            maxIterations, [], [], initialGuess);
        reconstruction.coefficients_tikhonov = coefficientsTik;
        reconstruction.sense_recon_stm_tikhonov = synthesize(coefficientsTik);
        reconstruction.sense_recon_stm_tik = reconstruction.sense_recon_stm_tikhonov;
        reconstruction.tikhonov = struct( ...
            'lambda', lambda, ...
            'flag', flagTik, ...
            'relativeResidual', relativeResidualTik, ...
            'iterations', iterationsTik, ...
            'residualHistory', residualHistoryTik, ...
            'elapsedTime', toc(tTik));
    end

    reconstruction.imageSize = coefficientSize(1:2);
    reconstruction.nFrame = op.nFrame;
    reconstruction.nBasis = op.nBasis;
    reconstruction.radialData = data;
    reconstruction.options = opts;

    % NRMSE calculations are intentionally retained for future ground truth.
    % if isfield(opts, 'groundTruth') && ~isempty(opts.groundTruth)
    %     reconstruction.NRMSE_stm = norm(opts.groundTruth(:) - ...
    %         reconstruction.sense_recon_stm(:)) / norm(opts.groundTruth(:));
    %     if logical(runTikhonov)
    %         reconstruction.NRMSE_stm_tikhonov = norm(opts.groundTruth(:) - ...
    %             reconstruction.sense_recon_stm_tikhonov(:)) / ...
    %             norm(opts.groundTruth(:));
    %     end
    % end

    function output = makeOutput(x, status, residual, count, history, elapsed, method)
        output = struct();
        output.coefficients = x;
        output.sense_recon_stm = synthesize(x);
        output.(method) = struct('flag', status, ...
            'relativeResidual', residual, 'iterations', count, ...
            'residualHistory', history, 'elapsedTime', elapsed);
    end

    function image = synthesize(x)
        image = reshape(op.synthesis(reshape(x, coefficientSize)), ...
            [coefficientSize(1), coefficientSize(2), op.nFrame]);
    end

    function value = getOption(options, name, defaultValue)
        if isfield(options, name) && ~isempty(options.(name))
            value = options.(name);
        else
            value = defaultValue;
        end
    end
end