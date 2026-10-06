function results = testReconstructGARSOSRadialZSlice()
%testRECONSTRUCTGARSOSRADIALZSLICE Test dimensions and centered z resizing.

    rng(17);
    nReadout = 4;
    nPartition = 5;
    nSpokes = 3;
    nCoil = 2;
    imageSizeZ = 8;
    gridSizeZ = 12;
    kernelSizeZ = 3;

    sampledData = complex(randn(nReadout*nPartition*nSpokes, nCoil), ...
        randn(nReadout*nPartition*nSpokes, nCoil));
    coordinates4 = zeros(nReadout, nPartition, nSpokes, 3);
    coordinates4(:, :, :, 1) = 0.1;
    coordinates4(:, :, :, 2) = -0.2;
    kz = reshape(linspace(-0.3, 0.3, nPartition), [1 nPartition 1]);
    coordinates4(:, :, :, 3) = repmat(kz, nReadout, 1, nSpokes);
    coordinates = reshape(coordinates4, [], 3);

    [sliceData, diagnostics] = reconstructGARSOSRadialZSlice( ...
        sampledData, coordinates, nReadout, nPartition, nSpokes, ...
        eye(3), imageSizeZ, gridSizeZ, kernelSizeZ, 4, ...
        'ReturnAllZ', true);

    if ~isequal(size(sliceData), [nReadout nSpokes nCoil])
        error('testReconstructGARSOSRadialZSlice:SliceSize', ...
            'The selected z slice has the wrong dimensions.');
    end
    if ~isequal(size(diagnostics.allZData), ...
            [nReadout nSpokes nCoil imageSizeZ])
        error('testReconstructGARSOSRadialZSlice:AllZSize', ...
            'The all-z output has the wrong dimensions.');
    end
    if diagnostics.nPartition ~= nPartition || ...
            diagnostics.gridSizeZ ~= gridSizeZ
        error('testReconstructGARSOSRadialZSlice:Diagnostics', ...
            'The reported z dimensions are inconsistent.');
    end

    resizeError = 0;
    for sizes = [12 8; 11 6; 8 11].'
        inputSize = sizes(1);
        outputSize = sizes(2);
        resizeInput = complex(randn(inputSize, 3), randn(inputSize, 3));
        expectedResize = shiftAndResize(resizeInput, [], ...
            [outputSize 3], 'ee', 'me');
        sourceIndices = shiftAndResize((1:inputSize).', [], ...
            outputSize, 'e', 'm');
        actualResize = zeros(outputSize, 3, 'like', resizeInput);
        valid = sourceIndices > 0;
        actualResize(valid, :) = resizeInput(sourceIndices(valid), :);
        thisError = norm(actualResize(:) - expectedResize(:)) / ...
            max(norm(expectedResize(:)), eps);
        resizeError = max(resizeError, thisError);
    end
    if resizeError > 1e-14
        error('testReconstructGARSOSRadialZSlice:Resize', ...
            'The centered z resize does not match shiftAndResize.');
    end

    results = struct('sliceSize', size(sliceData), ...
        'allZSize', size(diagnostics.allZData), ...
        'resizeRelativeError', resizeError);
end
