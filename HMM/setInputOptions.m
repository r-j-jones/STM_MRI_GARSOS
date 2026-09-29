function inputOptions = setInputOptions( ...
    frameParams, griddingParams, imageSize, gridSize, CA, Cb, visibleObjectSize)
%localBuildInputOptions Store variables in structs

    inputOptions = struct();
    inputOptions.imageSize = imageSize;
    inputOptions.gridSize = gridSize;
    inputOptions.nSpokesPerFrame = frameParams.nSpokesPerFrame;
    inputOptions.zSlice = frameParams.zSlice;
    inputOptions.calibrationSize = frameParams.calibrationSize;
    inputOptions.CoordinateTransform = CA;
    inputOptions.TranslationVector = Cb;
    inputOptions.VisibleObjectSize = visibleObjectSize;
    inputOptions.PhaseWidth = griddingParams.PhaseWidth;
    inputOptions.SpatialSigma = griddingParams.spatialSigma;
    inputOptions.DensityScale = griddingParams.DensityScale;
    inputOptions.UseGpu = griddingParams.UseGpu;
    inputOptions.DisplaySlice = griddingParams.DisplaySlice;
    inputOptions.DisplaySliceIndex = griddingParams.DisplaySliceIndex;
    inputOptions.ReturnDiagnostics = griddingParams.ReturnDiagnostics;
    inputOptions.KernelSize = griddingParams.kernelSize;
    inputOptions.OutputIsGpuArray = griddingParams.OutputIsGpuArray;
end