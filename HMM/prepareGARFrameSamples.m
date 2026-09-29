function outputs = prepareGARFrameSamples( ...
    sampleCoordinates, spokeIndex, nSpokesPerFrame, varargin)
%PREPAREGARFRAMESAMPLES Group GAR trajectory samples into temporal frames.
%
% USAGE:
% 
% outputs = prepareGARFrameSamples( ...
%     sampleCoordinates, spokeIndex, nSpokesPerFrame, ...
%     'CoordinateTransform', coordinateTransform);
% 
% Inputs:
%   sampleCoordinates   [Ns x 3] normalized [kx ky kz] coordinates.
%   spokeIndex          [Ns x 1] zero- or one-based spoke identifier.
%   nSpokesPerFrame     Positive integer number of consecutive spokes/frame.
%
% Name/value options:
%   'CoordinateTransform'
%       [3 x 3] trajectory transform. Default: eye(3).
%
% Output:
%   outputs             Scalar structure with four fields:
%
%       outputs.sampleCoordinates
%           1-by-nFrames cell array. Each cell contains the original
%           sample coordinates for that frame.
%
%       outputs.sampleMask
%           1-by-nFrames cell array. Each cell contains a logical Ns-by-1
%           mask selecting the samples belonging to that frame.
%
%       outputs.frameSampleIndices
%           1-by-nFrames cell array. Each cell contains the linear indices
%           of the samples belonging to that frame.
%
%       outputs.frameTrajectory
%           1-by-nFrames cell array. Each cell contains the transformed
%           trajectory coordinates for that frame.
%
% Only complete frames are returned. Any trailing spokes that do not form
% a complete frame are ignored.

    p = inputParser;

    addParameter(p, 'CoordinateTransform', eye(3), ...
        @(x) isnumeric(x) && isequal(size(x), [3 3]) && ...
        all(isfinite(x(:))));

    parse(p, varargin{:});

    validateattributes(sampleCoordinates, {'numeric'}, ...
        {'2d', 'ncols', 3, 'nonempty', 'finite'}, ...
        mfilename, 'sampleCoordinates');

    validateattributes(spokeIndex, {'numeric'}, ...
        {'vector', 'nonempty', 'finite'}, ...
        mfilename, 'spokeIndex');

    validateattributes(nSpokesPerFrame, {'numeric'}, ...
        {'scalar', 'integer', 'positive', 'finite'}, ...
        mfilename, 'nSpokesPerFrame');

    nSamples = size(sampleCoordinates, 1);
    spokeIndex = spokeIndex(:);

    if numel(spokeIndex) ~= nSamples
        error('prepareGARFrameSamples:SampleCountMismatch', ...
            ['sampleCoordinates and spokeIndex must contain the same ' ...
             'number of samples.']);
    end

    coordinateTransform = p.Results.CoordinateTransform;

    % Preserve the order in which spoke identifiers first occur.
    spokeValues = unique(spokeIndex, 'stable');

    % Include complete frames only.
    nFrames = floor(numel(spokeValues) / nSpokesPerFrame);

    if nFrames < 1
        error('prepareGARFrameSamples:InsufficientSpokes', ...
            ['There are not enough spokes to form one complete frame ' ...
             'containing %d spokes.'], nSpokesPerFrame);
    end

    spokeValues = spokeValues(1:nFrames * nSpokesPerFrame);

    outputs = struct();
    outputs.sampleCoordinates = cell(1, nFrames);
    outputs.sampleMask = cell(1, nFrames);
    outputs.frameSampleIndices = cell(1, nFrames);
    outputs.frameTrajectory = cell(1, nFrames);

    for iFrame = 1:nFrames
        frameSpokeRange = ...
            (iFrame - 1) * nSpokesPerFrame + (1:nSpokesPerFrame);

        frameSpokes = spokeValues(frameSpokeRange);

        sampleMask = ismember(spokeIndex, frameSpokes);
        frameSampleIndices = find(sampleMask);
        frameCoordinates = sampleCoordinates(sampleMask, :);
        frameTrajectory = frameCoordinates * coordinateTransform;

        outputs.sampleCoordinates{iFrame} = frameCoordinates;
        outputs.sampleMask{iFrame} = sampleMask;
        outputs.frameSampleIndices{iFrame} = frameSampleIndices;
        outputs.frameTrajectory{iFrame} = frameTrajectory;
    end
end