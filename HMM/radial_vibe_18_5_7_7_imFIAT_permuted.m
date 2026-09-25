function [manifest, manifestFile] = radial_vibe_18_5_7_7_imFIAT_permuted( ...
    metaManifest, reconstructionDir, referenceVolumeA, forceUpdate, savePermutedData)
% radial_vibe_18_5_7_7_imFIAT_permuted: Perform backprojection with gridding
% scanDir: the input directory
% metaManifest: metadata manifest
% reconstructionDir: the reconstruction directory
% referenceVolumeA: reference volume
% forceUpdate: boolean flag to force update

if nargin<5
    savePermutedData = false;
end
fprintf(' RUNNING: radial_vibe_18_5_7_7_imFIAT_permuted()-\n');
fprintf(' savePermutedData = %g\n',savePermutedData);

% Initialize timer
initTic = tic;

dstRootDir = reconstructionDir;
dstSeriesDirName = sprintf('%s_COMPLEX_PROJ_norm4', metaManifest.prefixString);
dstDir = fullfile(dstRootDir, dstSeriesDirName);
manifestFile = fullfile(dstDir, 'deid-manifest.mat');
reconstructionTimestamp = '2018-04-18 13:58:25 -04:00';

[manifestLoaded, manifest, ~] = loadIfExistAndTimestampNewerThan(manifestFile, reconstructionTimestamp);

if ~manifestLoaded || forceUpdate
    permissiveMakeDirIfNotExist(dstDir);    
    nSequences = numel(metaManifest.sequences);
    manifest = metaManifest;
   
    for iSequence = 1:nSequences

        seqTic = tic;

        fprintf('\n\n -- readVD11MultiRaidFileStructure, mapVBVD -- \n\n');
        s = readVD11MultiRaidFileStructure(metaManifest.sequences(iSequence).rawDataFile, 'CalibrationOnly', true);
        q = mapVBVD(metaManifest.sequences(iSequence).rawDataFile);
        if iscell(q)
            q = q{end};
        end
       
        fprintf('\n\n -- estimateCoilSensitivitieMaps2 -- \n\n');
        [noiseWhiteningTransform, channelSensitivityMaps, channelIds, ~] = estimateCoilSensitivitieMaps2(s);

        fprintf('\n\n -- prepareRadialVibeDataForGridding2_3 -- \n\n');
        [sampledData, kSpaceLocations, densityCompensation, iiSpokes, timeStamps, nPartitions, nSpokes, nSamples] = ...
            prepareRadialVibeDataForGridding2_3(q, s, noiseWhiteningTransform, channelIds);

        % [sampledData, kSpaceLocations, densityCompensation, iiSpokes, timeStamps, nPartitions, nSpokes, nSamples] = ...
        %     prepareRadialVibeDataForGridding2_3_1(q, s, noiseWhiteningTransform, channelIds);


       



        sensitivityMapArraySize = size(channelSensitivityMaps.A);
        periodicChannelSensitivityMaps = centeredArrayVolumeResize(channelSensitivityMaps, sensitivityMapArraySize(1,1:3).*[3 1 3]);
        periodicChannelSensitivityMaps.A = repmat(channelSensitivityMaps.A, [3 1 3 1]);
       
        nDimensions = 3;
        referenceSize = size(referenceVolumeA.A);
        if referenceSize(3) == 64
            imageSize = [224 224 floor(referenceSize(3) * 1.5)];
        else
            imageSize = [224 224 96];
        end
       
        gridSize = ceil(1.375 * imageSize);
        kernelSize = 7 * ones(1, nDimensions);
        nChannels = size(sampledData, 2);
        dicomImageSize = referenceSize;
        visibleObjectSize = referenceSize;
       
        L1 = [eye(3,3), floor(dicomImageSize'/2); 0 0 0 1];
        L2 = [diag([1 -1 -1]), [0; dicomImageSize(2) - 1; dicomImageSize(3) - 1]; 0 0 0 1];
        L3 = [referenceVolumeA.R * diag(referenceVolumeA.v), referenceVolumeA.r0; 0 0 0 1];
       
        Cs = L1 \ (L2 \ (L3 \ (eye(4,4) * L3 * L1)));
        CA = Cs(1:3,1:3);
        Cb = Cs(1:3,4);
       
        spatialSigma = 0; % px
        r0_machine = [s.m{2}.param.Meas{1}.DICOM.lGlobalTablePosSag, s.m{2}.param.Meas{1}.DICOM.lGlobalTablePosCor, s.m{2}.param.Meas{1}.DICOM.lGlobalTablePosTra];
       
        channelSensitivityMaps2 = periodicChannelSensitivityMaps;
        channelSensitivityMaps2.r0 = channelSensitivityMaps2.r0 + r0_machine';
       
        referenceVolumeC = centeredArrayVolumeResize(referenceVolumeA, imageSize);
        channelSensitivityMaps4 = resampleArrayVolumeToSame(channelSensitivityMaps2, referenceVolumeC, 'linear', 0);
        paddedChannelSensitivityArray = channelSensitivityMaps4.A;
       
        iiPhases = 0:1:nSpokes - 1;
        timeSeries = struct();
        timeSeries.iiSpokes = iiPhases;
        timeSeries.complexFile = fullfile(dstDir, sprintf('Sequence_%06.0f_ComplexData.raw', metaManifest.sequences(iSequence).sequenceId));
        timeSeries.imageSize = imageSize;
        timeSeries.dicomImageSize = dicomImageSize;
        timeSeries.nSpokes = nSpokes;
        timeSeries.sensitivityMapArraySize = sensitivityMapArraySize;
        % timeSeries.kSpaceLocations = kSpaceLocations;
        % timeSeries.dc = densityCompensation;


        %% Check density comp values

        % check_dcf( densityCompensation, nPartitions, nSpokes, nSamples );


        %% STM DATA PREP

        % Inputs:
        %   sampledData          [Ns x Nc] complex corrected radial samples.
        %   sampleCoordinates    [Ns x 3] normalized [kx ky kz] coordinates.
        %   densityCompensation  [Ns x 1] radial density-compensation weights.
        %   spokeIndex           [Ns x 1] zero- or one-based spoke identifier.
        %   sensitivityMaps      [Nx x Ny x Nz x Nc] complex maps, or an object with
        %                        an A property containing that array.
        %   imageSize            [1 x 3] reconstructed volume dimensions.
        %   gridSize             [1 x 3] oversampled Cartesian grid dimensions.
        %   nSpokesPerFrame      positive integer number of consecutive spokes/frame.
        %   zSlice               one-based image-space z index.
        %   calibrationSize      [NxCal NyCal] central rectangular k-space size.
        %
        % Name/value options:
        %   'KernelSize'         [1 x 3], default [7 7 7].
        %   'OutputIsGpuArray'   logical, default false.
        %   'CoordinateTransform' [3 x 3] HMM trajectory transform.
        %   'TranslationVector'  [3 x 1] HMM translation vector.
        %   'VisibleObjectSize'  [1 x 3] HMM visible object size.
        %   'PhaseWidth'         scalar HMM phase-width normalization.
        %   'SpatialSigma'       scalar HMM Gaussian weighting parameter.
        %   'DensityScale'       additional scalar after HMM density scaling.
        %   'ReturnDiagnostics'  logical, default false.


        nSpokesPerFrame = 6;
        zSlice = 1:imageSize(3);
        calibrationSize = [32 32];


        [kCal, outputs, diagnostics] = reconstructGARFrameForSTM( ...
            sampledData, kSpaceLocations, densityCompensation, iiSpokes, ...
            channelSensitivityMaps4, imageSize, gridSize, nSpokesPerFrame, zSlice, ...
            calibrationSize, ...
            'CoordinateTransform', CA, ...
            'TranslationVector', Cb, ...
            'VisibleObjectSize', visibleObjectSize, ...
            'PhaseWidth', 1, ...
            'SpatialSigma', spatialSigma, ...
            'DensityScale', 1, ...
            'UseGpu', true, ...
            'DisplaySlice', false, ...
            'ReturnDiagnostics', true, ...
            'KernelSize', kernelSize);



        %% ORIGINAL HMM CODE BELOW:
        if 0 == 1
       
            nNominalSpokesInBlock = 24;        
            nBlocks = round(numel(iiPhases) / nNominalSpokesInBlock);
            jjSpokesInBlocks = partitionAllTasks(1:numel(iiPhases), nBlocks);
   
            griddingParams = struct();
            griddingParams.nDimensions = nDimensions;
            griddingParams.referenceSize = referenceSize;
            griddingParams.imageSize = imageSize;
            griddingParams.gridSize = gridSize;
            griddingParams.kernelSize = kernelSize;
            griddingParams.nChannels = nChannels;
            griddingParams.dicomImageSize = dicomImageSize;
            griddingParams.visibleObjectSize = visibleObjectSize;
            griddingParams.spatialSigma = spatialSigma;
            griddingParams.r0_machine = r0_machine;
            griddingParams.nNominalSpokesInBlock = nNominalSpokesInBlock;
           
            switch_gpu(1, 2);
           
            fprintf('\n\n -- Storing data in cell arrays -- \n\n');
            sampledDatas = cell(nBlocks, 1);
            kSpaceLocationss = cell(nBlocks, 1);
            densityCompensations = cell(nBlocks, 1);
           
            for iBlock = 1:nBlocks
                jjSpokesInBlock = jjSpokesInBlocks{iBlock};
                iiSpokesInBlock = iiPhases(jjSpokesInBlock);
                [mKeep, ~] = ismember(iiSpokes, iiSpokesInBlock);
                sampledDatas{iBlock} = sampledData(mKeep, :);
                kSpaceLocationss{iBlock} = kSpaceLocations(mKeep, :);
                densityCompensations{iBlock} = densityCompensation(mKeep, :);
            end
           
            outputFilePath = timeSeries.complexFile;
            clear sampledData kSpaceLocations densityCompensation
            % outputLogFile = fullfile(dstDir, 'fwrite_bit_log.txt');
            % idlog = fopen(outputLogFile,'w');
   
   
            fprintf('\n\n -- Starting parfor recon -- \n\n');
           
            parfor iBlock = 1:nBlocks
                tBlock = tic;
                phaseWidth = 1;
                jjSpokesInBlock = jjSpokesInBlocks{iBlock};
                nSpokesInBlock = numel(jjSpokesInBlock);
               
                s0 = sampledDatas{iBlock};
                k0 = tprod(kSpaceLocationss{iBlock}, [1 -5], CA, [-5 2]);
                c0 = exp(1i * 2 * pi * tprod(kSpaceLocationss{iBlock}, [1 -5], Cb, [-5 2]));
               
                d20 = sum(k0.^2, 2);
                w20 = exp(-2 * pi * (d20 * spatialSigma.^2));
                w0 = densityCompensations{iBlock} * prod(visibleObjectSize(1:2)) * 2;
                w0f = w0 / phaseWidth;
                sw0 = bsxfun(@times, s0, w0f .* w20 .* c0);
               
                k03 = mat2cell(k0, size(k0, 1) / nSpokesInBlock * ones(1, nSpokesInBlock), nDimensions);
                sw03 = mat2cell(sw0, size(k0, 1) / nSpokesInBlock * ones(1, nSpokesInBlock), nChannels);
                S03pre = cell(size(sw03));
                paddedChannelSensitivityArray03 = gpuArray(paddedChannelSensitivityArray);
               
                for i03Spoke = 1:nSpokesInBlock
                    G = GpuReconstructur20181008(gridSize, imageSize, kernelSize, gpuArray(k03{i03Spoke}), paddedChannelSensitivityArray03);
                    S03pre{i03Spoke} = G.grid(sw03{i03Spoke}, false);
                end
               
                I0 = cat(4, S03pre{:});
                if savePermutedData
                    I0 = permute(I0, [1 2 4 3]);
                end
   
                fidComplex = fopen(outputFilePath, 'r+');            
                fseek(fidComplex, prod(imageSize) * 4 * (jjSpokesInBlock(1) - 1), -1);
                fwrite(fidComplex, real(I0), 'single');
                fseek(fidComplex, prod(imageSize) * 4 * (nSpokes + jjSpokesInBlock(1) - 1), -1);
                fwrite(fidComplex, imag(I0), 'single');        
                fclose(fidComplex);
   
                % fprintf(idlog,'spokes %d %d, real_start %d, imag_start %d\n',...
                %     jjSpokesInBlock(1), jjSpokesInBlock(end), ...
                %     prod(imageSize) * 4 * (jjSpokesInBlock(1) - 1), ...
                %     prod(imageSize) * 4 * (nSpokes + jjSpokesInBlock(1) - 1));
   
   
               
                fprintf('Reconstructed spokes %.0f through %.0f in %.2f seconds.\n', jjSpokesInBlock(1), jjSpokesInBlock(end), toc(tBlock));
                fprintf('Current time: %s\n', iso8601Now());
            end
   
            % fclose(idlog);
   
   
            seqToc = toc(seqTic);
           
            timeSeries.dataIsPermuted = savePermutedData;
            timeSeries.griddingParams = griddingParams;
            manifest.sequences(iSequence).timeSeries = timeSeries;
            manifest.sequences(iSequence).dataIsPermuted = savePermutedData;
            manifest.sequences(iSequence).elapProcessingTime = seqToc;
        end
       
        for iSequence = 1:nSequences
            assert(exist(manifest.sequences(iSequence).timeSeries.complexFile, 'file') == 2);
        end
   
        elapTime = toc(initTic);
        manifest.totalElapsedTime = elapTime;
       
        % Display total elapsed time
        fprintf('Total elapsed time: %.2f seconds.\n', elapTime);
       
        saveWithTimestamp(manifestFile, manifest, iso8601Now());
        [manifestLoaded, manifest, ~] = loadIfExistAndTimestampNewerThan(manifestFile, reconstructionTimestamp);
        assert(manifestLoaded);

    end
   
end



end