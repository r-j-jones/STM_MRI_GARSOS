function plans = buildGARSOSRadialNUFFTPlans(kxy, imageSize, varargin)
%buildGARSOSRADIALNUFFTPLANS Build one pure 2-D NUFFT plan per frame.
%
% kxy: [nReadout, nSpokes, nFrame], normalized kx+i*ky coordinates.
% The underlying Fessler NUFFT receives radians; this function performs the
% required multiplication by 2*pi exactly once.

    p = inputParser;
    addParameter(p, 'Jd', [6 6], @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'Kd', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 2));
    parse(p, varargin{:});

    validateattributes(kxy, {'numeric'}, {'3d', 'nonempty'});
    validateattributes(imageSize, {'numeric'}, {'vector', 'numel', 2, 'integer', 'positive'});
    if any(~isfinite(kxy(:)))
        error('buildGARSOSRadialNUFFTPlans:NonFiniteTrajectory', ...
            'The normalized trajectory contains nonfinite values.');
    end

    nFrame = size(kxy, 3);
    if isempty(p.Results.Kd)
        Kd = floor(1.5 * imageSize);
    else
        Kd = p.Results.Kd(:).';
    end
    Nd = imageSize(:).';
    Jd = p.Results.Jd(:).';
    plans = repmat(struct('st', [], 'nReadout', [], 'nSpokes', [], ...
        'normalizedTrajectory', [], 'omega', [], 'imageSize', Nd, ...
        'trajectoryUnits', 'cycles/FOV-like normalized units', ...
        'omegaUnits', 'radians'), 1, nFrame);

    for iFrame = 1:nFrame
        kFrame = kxy(:, :, iFrame);
        omega = [real(kFrame(:)), imag(kFrame(:))] * 2 * pi;
        plans(iFrame).st = nufft_init( ...
            omega, Nd, Jd, Kd, Nd / 2, 'kaiser');
        plans(iFrame).nReadout = size(kxy, 1);
        plans(iFrame).nSpokes = size(kxy, 2);
        plans(iFrame).normalizedTrajectory = kFrame;
        plans(iFrame).omega = omega;
    end
end
