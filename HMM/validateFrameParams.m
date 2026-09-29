function frameParams = validateFrameParams(frameParams)
%validateFrameParams Validate frameParams

    if ~isstruct(frameParams) || ~isscalar(frameParams)
        error('frameParams must be a scalar struct.');
    end
    required = {'nSpokesPerFrame', 'zSlice', 'calibrationSize'};
    if ~all(isfield(frameParams, required))
        error('frameParams must contain nSpokesPerFrame, zSlice, calibrationSize.');
    end
    validateattributes(frameParams.nSpokesPerFrame, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    if ~isempty(frameParams.zSlice)
        validateattributes(frameParams.zSlice, {'numeric'}, ...
            {'vector', 'integer', 'positive'});
    end
    validateattributes(frameParams.calibrationSize, {'numeric'}, ...
        {'vector', 'numel', 2, 'integer', 'positive'});
end