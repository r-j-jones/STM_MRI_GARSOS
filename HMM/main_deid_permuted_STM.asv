% main_deid_permuted_STM(inputPath, slow_state_num, breath_state_num, gi_state_num, HTres_sigmas, train_splits, opts)
% 
% RJ  |  09-20-2026
% 
% This runs MMM for KRR manuscript:
%   - specify train/test split to use
%   - specify HTres sigmas to use
%   - will use new "permute" functions for faster .raw reading/writing
%   - will use updated deformable alignment
%   
% EXAMPLE USAGE:
% 
% inputPath = '/mnt/ibrixfs04-Kspace/motion_patients/P0004/E01';
% breath_state_num = 21;
% HTres_sigmas = [
%     1,1;
%     1,3;
%     ];
% train_splits = {'1-1000','1001-2000'};
% opts.forceRaw = false;
% opts.forceUpdate = false;




% inputPath = '/mnt/ibrixfs04-Kspace/motion_patients/P0021/E02';
% outputPath = '/RadOnc-MRI1/Student_Folder/rjones/STM/InitialTest_09-20-2026';
% opts = [];
% 
% cd('/RadOnc-MRI1/Student_Folder/rjones/STM/STM_MRI_GARSOS/HMM');


function main_deid_permuted_STM(inputPath, outputPath, opts)

    fprintf('\n-running main_deid_permuted()-\n  [ %s ]  \n\n%s\n\n',datetime,inputPath);

    if ~isdeployed
        addpath('/RadOnc-MRI1/Student_Folder/rjones/STM/STM_MRI_GARSOS');
        addpath('/RadOnc-MRI1/Student_Folder/rjones/STM/STM_MRI_GARSOS/HMM');
        % addpath('/RadOnc-MRI1/Student_Folder/rjones/MR_Motion_Modeling/new_source_code/YZ_data');
        % addpath('/RadOnc-MRI1/Student_Folder/rjones/KernelRegression/code/version_04-15-2024/main-utils');
        % %%% Set paths
        % set_main_MMM_paths();
    end

    if nargin<3
        opts.forceRaw = false;
        opts.forceUpdate = false;
        opts.runPostBreathing = false;
        opts.viz = true;
        opts.paropts.parallelWrite = false;
        opts.paropts.nParforWorkers = 12;
    else
        if ~isfield(opts,'forceRaw'), opts.forceRaw = false; end
        if ~isfield(opts,'forceUpdate'), opts.forceUpdate = false; end
        if ~isfield(opts,'viz'), opts.viz = true; end
        if ~isfield(opts,'paropts')
            opts.paropts.parallelWrite = false; 
            opts.paropts.nParforWorkers = 12; 
        end
    end
    
    fprintf('\n\n----------------------------------')
    fprintf('\n Running main_deid_permuted_STM \n');
    fprintf('\n----------------------------------\n\n');

    viz = opts.viz;
    paropts = opts.paropts;
    forceUpdate = opts.forceUpdate;

    %[ Get metadata on current run  
    %hostname
    [~, hostname] = system('hostname');
    hostname = strtrim(hostname); % Remove trailing newline
    startDatetime = datetime;    
    %matlab-version
    [matlabVersion.v,matlabVersion.d] = version;

    %% (0) Setup + prep

    %[ Get subject deid nb + experiment nb
    [subjNbStr, experimentNbStr, scanStr] = parseDeidInputPath( inputPath );
    

    %[ Get input path
    reconstructionDir = fullfile(inputPath,'reconstruction');
    load(fullfile(reconstructionDir,'input.mat'),'kSpaceID1', 'kSpaceID2');
    scanDir = inputPath;
    rawDir = fullfile(inputPath, 'raw');

    %[ Get output path
    outputDir = fullfile(outputPath,'hmm-data');
    if ~exist(outputDir,'dir')
        mkdir(outputDir); 
    end


    % %[ create dirs
    % plotDir = fullfile(reconstructionDir,'plots-motion-model');
    % if ~exist(plotDir,'dir'), mkdir(plotDir); end
    % logFileDir = fullfile(reconstructionDir,'model-log-files');
    % if ~exist(logFileDir,'dir'), mkdir(logFileDir); end
    % globalPlotDir = '/RadOnc-MRI1/Student_Folder/rjones/MR_Motion_Modeling/motion-deid-plots/BreathingMS-movies';
    % if ~exist(globalPlotDir,'dir'), mkdir(globalPlotDir); end

    %[ Save 
    fInputDeid = fullfile(outputDir,'deid-input.mat');
    if ~exist(fInputDeid,'file') || forceUpdate
        save(fInputDeid,'inputPath',...
            'reconstructionDir','opts','scanDir','hostname','startDatetime',...
            'subjNbStr','experimentNbStr','scanStr','rawDir','matlabVersion');
    end


    %%% Load reference Image here
    if isempty(gcp('nocreate')), parpool(12); end
    [referenceImage, templateFile, ~] = loadReferenceDicomImage( scanDir, true );
    % if viz 
    %     if ~exist(fullfile(plotDir,'im_referenceImage.png'),'file')
    %         displayReferenceArrayVolume( referenceImage, plotDir ); 
    %     else
    %         displayReferenceArrayVolume( referenceImage );
    %     end
    % end

    %%% K-space ID identification
    [kSpaceID, ~] = parseKspaceIds( kSpaceID1, kSpaceID2 );    

    %%% To store manifest file paths
    fmanifests = []; fcheckpoints = [];


    %% (1) Load k space data

    [metaManifest, fmanifests.meta] = radial_vibe_18_5_4_timestamp_imFIAT( ...
        kSpaceID, scanDir, outputDir, opts.forceUpdate ); % Meta

    %% Get acquisition parameters struct

    [acqparamsManifest, fmanifests.acqparams] = radial_vibe_get_acq_params( ...
        kSpaceID, scanDir, outputDir, opts.forceUpdate);
    
    %% (2) Back project and regrid k space data to image space

    [projManifest, fmanifests.proj] = radial_vibe_18_5_7_7_imFIAT_permuted( ...
        metaManifest, outputDir, referenceImage, opts.forceUpdate, false);


    disp('done getting stm data...');



    %% (3) 3d distortion correction

    [dis3dProjManifest, fmanifests.dis3d_proj] = radial_vibe_28_dis3d_proj_2_imFIAT_permuted( ...
        projManifest, reconstructionDir, referenceImage, templateFile, opts.forceUpdate, true, 25);

    % fplotDisco3d = fullfile(plotDir,'dis3d_back-proj_subplots.png');
    % savePlotFlag = false; 
    % if ~exist(fplotDisco3d,'file') || forceUpdate, savePlotFlag=true; end
    % [~,~] = display_dis3d_backproj_slice( dis3dProjManifest, 1:25, 50, plotDir, savePlotFlag ); 
    
    % %[ checkpoint - Save the timers, manifests
    % fcheckpoints.log1 = fullfile(logFileDir,'main_deid_perumted-checkpoint_post-disco3d.mat');
    % if ~exist(fcheckpoints.log1,'file') || forceUpdate
    %     save(fcheckpoints.log1, 'fmanifests', 'projManifest', 'dis3dProjManifest', 'metaManifest');
    % end
    
    %% (4) Reconstruct high temporal resolution images, still use padding for this one
    
    HTres_stride = 1;
    HTres_params.stride = HTres_stride;
    HTres_params.HTres_sigmas = HTres_sigmas;
    HTmanifests = [];
    % HTtimers = [];
    for ht_ind = 1:size(HTres_sigmas,1)
        HTsigma_min = HTres_sigmas(ht_ind,1);
        HTsigma_max = HTres_sigmas(ht_ind,2);
        fprintf('smin = %d, smax = %d\n',HTsigma_min, HTsigma_max);
        sigmasStr = sprintf('sigma_%d_%d',HTsigma_min,HTsigma_max);
        
        % atic = tic; acpu = cputime;

        [torandoManifestHT, fManifestHT] = radial_vibe_28_combine_spoke_volumes_imFIAT_permuted( ...
            dis3dProjManifest, reconstructionDir, HTsigma_min, HTsigma_max, HTres_stride, opts.forceUpdate, ...
            'HT-new', referenceImage, paropts); 

        % atoc = toc(atic); acpu = cputime - acpu;

        HTmanifests.(sigmasStr) = torandoManifestHT;
        % HTtimers.(sigmasStr).elap = atoc; %torandoManifestHT.totalProcessingTime;
        % HTtimers.(sigmasStr).cpu = acpu;
        fmanifests.ht_rec.(sigmasStr) = fManifestHT;
    end
    

    % %[ checkpoint - Save the timers, manifests
    % fcheckpoints.log2 = fullfile(logFileDir,'main_deid_perumted-checkpoint_post-HT-recon.mat');
    % if ~exist(fcheckpoints.log2,'file') || forceUpdate
    %     save(fcheckpoints.log2, 'fmanifests', 'projManifest', 'dis3dProjManifest', 'metaManifest',...
    %         'HTmanifests','HTtimers','HTres_params');
    % end


    %% (5) Rigid registration of images    

    % Loop through HT-res sigmas, load/perform rigid reg
    htregTransforms = [];
    % htregTimers = [];
    for ht_ind = 1:size(HTres_sigmas,1)
        HTsigma_min = HTres_sigmas(ht_ind,1);
        HTsigma_max = HTres_sigmas(ht_ind,2);
        fprintf('smin = %d, smax = %d\n',HTsigma_min, HTsigma_max);
        sigmasStr = sprintf('sigma_%d_%d',HTsigma_min,HTsigma_max);

        % atic = tic; acpu = cputime;
        [registrationTransform, fRegistrationTransform] = radial_vibe_11_register_21_spokes_2_6_imFIAT(...
            scanDir, HTmanifests.(sigmasStr), reconstructionDir, opts.forceUpdate, referenceImage, paropts); %s_REG    
        % atoc = toc(atic); acpu = cputime - acpu;
        
        % htregTimers.(sigmasStr).elap = atoc;
        % htregTimers.(sigmasStr).cpu = acpu;

        htregTransforms.(sigmasStr) = registrationTransform;
        fmanifests.ht_reg.(sigmasStr) = fRegistrationTransform;

        % %[ Compare transforms from 1st and 2nd rounds of rigid registration
        % evalInfo = registrationTransform.evaluationInformation;
        % [tfmDiffs, SI_tfmDiffs] = compareEvaluationInformation(evalInfo{2}, evalInfo{3});
        % htregTransforms.(sigmaStr).tfmDiffs = tfmDiffs;
        % htregTransforms.(sigmaStr).SI_tfmDiffs = SI_tfmDiffs;

    end
    
    % %[ checkpoint - Save the timers, manifests
    % fcheckpoints.log3 = fullfile(logFileDir,'main_deid_perumted-checkpoint_post-HT-rigid-reg.mat');
    % if ~exist(fcheckpoints.log3,'file') || forceUpdate
    %     save(fcheckpoints.log3, 'fmanifests', 'projManifest', 'dis3dProjManifest', 'metaManifest',...
    %         'HTmanifests','HTtimers','HTres_params','htregTimers');
    % end

    %% (6) Extract breathing motion signal
    
    % Loop through HT-res sigmas, load/perform rigid reg
    zsiManifests = []; 
    % zsiTimers = [];
    htresSigmasCell = cell(1,size(HTres_sigmas,1));
    for ht_ind = 1:size(HTres_sigmas,1)
        HTsigma_min = HTres_sigmas(ht_ind,1);
        HTsigma_max = HTres_sigmas(ht_ind,2);
        fprintf('smin = %d, smax = %d\n',HTsigma_min, HTsigma_max);
        sigmasStr = sprintf('sigma_%d_%d',HTsigma_min,HTsigma_max);
        currTransforms = htregTransforms.(sigmasStr);
        htresSigmasCell{ht_ind} = sigmasStr;

        % atic = tic; acpu = cputime;

        [motionSignalManifest,fmotionSignalManifest] = radial_vibe_28_motion_signal_2_imFIAT(...
            scanDir, projManifest, metaManifest, currTransforms, reconstructionDir, ...
            opts.forceUpdate, referenceImage, sigmasStr);

        % atoc = toc(atic); acpu = cputime - acpu;

        % zsiTimers.(sigmasStr).elap = atoc;
        % zsiTimers.(sigmasStr).cpu = acpu;
        zsiManifests.(sigmasStr) = motionSignalManifest;
        fmanifests.zsi.(sigmasStr) = fmotionSignalManifest;
    end
    
    % if 0 == 1
    %     %[ Plot zSI breathing motion signal traces (1/3),(2/5),(orig 2/5)
    %     compare_zSI_traces( reconstructionDir, zsiManifests, plotDir); 
    %     compare_zSI_traces_general(reconstructionDir, zsiManifests, plotDir, htresSigmasCell);
    % end


    % %[ checkpoint - Save the timers, manifests
    % fcheckpoints.log4 = fullfile(logFileDir,'main_deid_perumted-checkpoint_post-zSI.mat');
    % if ~exist(fcheckpoints.log4,'file') || forceUpdate
    %     save(fcheckpoints.log4, 'fmanifests', 'projManifest', 'dis3dProjManifest', 'metaManifest',...
    %         'HTmanifests','HTtimers','HTres_params','htregTimers', 'zsiTimers', 'zsiManifests');
    % end

    %% (7) 21 B-MS; Canonical Breathing Motion State (B-MS) view-sharing (VS) image reconstruction 
    
    % %%% Train split: Canonical Breathing Motion State Tornado Image Recon
    % train_spoke_info = parseTrainSplitString(train_splits);
    % nSplits = length(train_splits);

    %[ Loop over {11,21} MS's reconed
    for nmsInd=1:nMotionStateNbs
        nbOfBms = nBreathMsToRecon(nmsInd);
        nbOfBms_tag = sprintf('ms%d',nbOfBms);
    
        % %[ Loop over train splits
        % for splitInd = 1:nSplits
        %     split_info = train_spoke_info{splitInd};
        %     split_inds = split_info.inds;
        %     split_tag = split_info.tag;
        %     split_str = sprintf('split%d',splitInd);
        %     fprintf('-[%s]-\n -Processing: %s\n',datetime,split_tag);
    
        %[ Loop over HT res recons
        for ht_ind = 1:size(HTres_sigmas,1)
            HTsigma_min = HTres_sigmas(ht_ind,1);
            HTsigma_max = HTres_sigmas(ht_ind,2);
            fprintf('smin = %d, smax = %d\n',HTsigma_min, HTsigma_max);
            sigmasStr = sprintf('sigma_%d_%d',HTsigma_min,HTsigma_max); 
            suffixStr = sprintf('_HT-%s',sigmasStr);

            %[ Get the current zSI brething motion signal
            zsiSignal = zsiManifests.(sigmasStr).zSI;
           
            %[ Set view-sharing Gaussian sigmas(min/max filter widths) + stride
            [filterParams,~] = compute_BMS_filter_params( 1:length(zsiSignal), nbOfBms, zsiSignal );
            % assert( all( zsiTrain(:)==zsiSignal(:)));

            % %[ do canonical breathing MS VS recon (& time)
            % % atic = tic; acpu = cputime;
            % % [breathingTimeSeries,fBMSmanifestFile] = radial_vibe_28_combine_sorted_spoke_volumes_imFIAT_Split(...
            % %     dis3dProjManifest, zsiTrain, reconstructionDir, filterParams.smin, filterParams.smax, ...
            % %     filterParams.stride, opts.forceUpdateBmsRec, suffixStr, referenceImage, split_inds, true);
            % % atoc = toc(atic); acpu = cputime - acpu;

            % if nmsInd==1 && ht_ind==1
            %     onlyWriteBmsManifest = true;
            % else
            %     onlyWriteBmsManifest = false;
            % end
            onlyWriteBmsManifest = false;

            atic = tic; acpu = cputime;
            [breathingTimeSeries,fBMSmanifestFile,~] = radial_vibe_28_combine_sorted_spoke_volumes_imFIAT_Full(...
              dis3dProjManifest, zsiSignal, reconstructionDir, filterParams.smin, filterParams.smax, ...
                filterParams.stride, opts.forceUpdateBmsRec, suffixStr, referenceImage, true, onlyWriteBmsManifest);
            atoc = toc(atic); acpu = cputime - acpu;

            %[ Store results
            bmsReconTimers.(nbOfBms_tag).(sigmasStr).elap = atoc;
            bmsReconTimers.(nbOfBms_tag).(sigmasStr).cpu = acpu;
            bmsReconManifests.(nbOfBms_tag).(sigmasStr) = breathingTimeSeries;
            fmanifests.bms_rec.(nbOfBms_tag).(sigmasStr) = fBMSmanifestFile;

            %[ Make gifs
            filedir = fileparts(breathingTimeSeries.niftiFiles{1});
            gifname = sprintf('BMS-movie_Full_%s%s_%s',scanStr,suffixStr,nbOfBms_tag);
            slicenbs=[];
            gifParams=[];
            make_deid_BMS_gifs( filedir, globalPlotDir, gifname, slicenbs, gifParams );

        end
       
    end

    fprintf('\n\n all done, now deformable reg next..\n\n');

    

    if 0 == 1
        % files = dir([reconstructionDir filesep 'Tornado5_MS_*']);
        % isDir = [files.isdir];                       % Logical index for directories
        % % Exclude '.' and '..'
        % names = {files(isDir).name};
        % names = names(~ismember(names, {'.', '..'}));
        % filtered = names(~contains(names, 'HT-sigma'));
        % filtered = filtered(~contains(filtered, 'DEF'));
    
        %[ Make gifs
        cmd = ['ls -d ' reconstructionDir '/Tornado5_MS_* | grep -v HT-sigma | grep -v DEF'];
        [~,fullBmsFiledir]=system(cmd,'-echo');
        fullBmsFiledir = strtrim(fullBmsFiledir);
        gifname = sprintf('BMS-movie_full_%s',scanStr);
        slicenbs=[];
        gifParams=[];
        make_deid_BMS_gifs( fullBmsFiledir, globalPlotDir, gifname, slicenbs, gifParams );
    end

    % if opts.onlySaveBmsGifs
    %     fprintf(' -only saving BMS gifs - exiting!\n\n');
    %     return
    % end


    % %[ checkpoint - Save the timers, manifests
    % fcheckpoints.log5 = fullfile(logFileDir,'main_deid_perumted-checkpoint_post-BMS-rec.mat');
    % if ~exist(fcheckpoints.log5,'file') || opts.forceUpdateBmsRec %forceUpdate
    %     save(fcheckpoints.log5, 'fmanifests', 'projManifest', 'dis3dProjManifest', 'metaManifest',...
    %         'HTmanifests','HTtimers','HTres_params','htregTimers', 'zsiTimers', 'zsiManifests',...
    %         'bmsReconTimers','bmsReconManifests');
    % end

    %% ((( 8 ))) Deformable registration of canonical breathing motion state volumes w/ niftyreg

    % % Deformable registration of breathing motion state images
    
    niftyregTimers = [];
    niftyregManifests = [];

    %[ Loop over {11,21} MS's reconed
    for nmsInd=1:nMotionStateNbs
        nbOfBms = nBreathMsToRecon(nmsInd);
        % nbOfBms_str = sprintf('bms_%d',nbOfBms);
        nbOfBms_tag = sprintf('ms%d',nbOfBms);

        bmsReconManifest = bmsReconManifests.(nbOfBms_tag);

        % %[ Loop over train splits
        % for splitInd = 1:nSplits
        %     split_info = train_spoke_info{splitInd};
        %     % split_inds = split_info.inds;
        %     split_tag = split_info.tag;
        %     split_str = sprintf('split%d',splitInd);
        %     fprintf('-[%s]-\n -Processing: %s\n',datetime,split_tag);
    
        %[ Loop over HT res recons
        for ht_ind = 1:size(HTres_sigmas,1)
            HTsigma_min = HTres_sigmas(ht_ind,1);
            HTsigma_max = HTres_sigmas(ht_ind,2);
            
            sigmasStr = sprintf('sigma_%d_%d',HTsigma_min,HTsigma_max); 
            % suffixStr = sprintf('_%s_HT-%s',split_tag,sigmasStr);

            %[ Get motion state VS recon manifest
            motionStateManifest = bmsReconManifest.(sigmasStr);

            fprintf('--- SUBMITTING NIFTYREG DIR OF BMS FOR: ---\n');
            fprintf(' Nb of BMS = %d\n',nbOfBms);
            % fprintf(' Train split = %s (%s)\n', split_tag, split_str);
            fprintf(' HT-smin = %d, HT-smax = %d\n',HTsigma_min, HTsigma_max);

            %[ do canonical breathing MS deformable reg (niftyreg) (& time)
            atic = tic; acpu = cputime;
            [~, ~, bmsRegManifest, bmsRegManifestFile] = ...
                runMotionStateDeformableRegistration_Split( scanDir, motionStateManifest, ...
                reconstructionDir, true, referenceImage, true, false);
            atoc = toc(atic); acpu = cputime - acpu;

            %[ Store results
            niftyregTimers.(nbOfBms_tag).(sigmasStr).elap = atoc;
            niftyregTimers.(nbOfBms_tag).(sigmasStr).cpu = acpu;
            niftyregManifests.(nbOfBms_tag).(sigmasStr) = bmsRegManifest;
            fmanifests.bms_reg.(nbOfBms_tag).(sigmasStr) = bmsRegManifestFile;

        end
           
        % end
    end

    % %[ checkpoint - Save the timers, manifests
    % fcheckpoints.log6 = fullfile(logFileDir,'main_deid_perumted-checkpoint_post-BMS-reg.mat');
    % if ~exist(fcheckpoints.log6,'file') || opts.forceUpdateBmsRec %forceUpdate
    %     save(fcheckpoints.log6, 'fmanifests', 'projManifest', 'dis3dProjManifest', 'metaManifest',...
    %         'HTmanifests','HTtimers','HTres_params','htregTimers', 'zsiTimers', 'zsiManifests',...
    %         'bmsReconTimers','bmsReconManifests','niftyregTimers','niftyregManifests');
    % end

    %% ( Last ) finish up

    finishDatetime = datetime;
    % if forceUpdate || exist(fInputDeid,'file') || opts.forceUpdateBmsRec 
    %     save(fInputDeid,'finishDatetime','-append');
    % end

    fprintf('\n\n-FINISHED running main_deid_permuted()-\n  [ %s ] \n [ %s ]  \n\n',scanStr,finishDatetime);




    % if 0 == 1
    %     %% [[[ 8: Analyze B-MS deformable registrations
    % 
    %     %[ Loop over {11,21} MS's reconed
    %     for nmsInd=1:nMotionStateNbs
    %         nbOfBms = nBreathMsToRecon(nmsInd);
    %         % nbOfBms_str = sprintf('bms_%d',nbOfBms);
    %         nbOfBms_tag = sprintf('ms%d',nbOfBms);
    % 
    %         bmsReconManifest = bmsReconManifests.(nbOfBms_tag);
    % 
    %         %[ Loop over train splits
    %         for splitInd = 1:nSplits
    %             split_info = train_spoke_info{splitInd};
    %             % split_inds = split_info.inds;
    %             split_tag = split_info.tag;
    %             split_str = sprintf('split%d',splitInd);
    %             % fprintf('-[%s]-\n -Processing: %s\n',datetime,split_tag);
    % 
    %             %[ Loop over HT res recons
    %             for ht_ind = 1:size(HTres_sigmas,1)
    %                 HTsigma_min = HTres_sigmas(ht_ind,1);
    %                 HTsigma_max = HTres_sigmas(ht_ind,2);
    % 
    %                 sigmasStr = sprintf('sigma_%d_%d',HTsigma_min,HTsigma_max); 
    %                 % suffixStr = sprintf('_%s_HT-%s',split_tag,sigmasStr);
    % 
    %                 %[ Get motion state VS recon manifest
    %                 motionStateManifest = bmsReconManifest.(split_str).(sigmasStr);
    % 
    %                 fprintf('--- ANALYZING NIFTYREG DIR OF BMS FOR: ---\n');
    %                 fprintf(' Nb of BMS = %d\n',nbOfBms);
    %                 fprintf(' Train split = %s (%s)\n', split_tag, split_str);
    %                 fprintf(' HT-smin = %d, HT-smax = %d\n',HTsigma_min, HTsigma_max);
    % 
    %                 %[ do canonical breathing MS deformable reg (niftyreg) (& time)
    %                 [~, ~, bmsRegManifest, ~] = ...
    %                     runMotionStateDeformableRegistration_Split( scanDir, motionStateManifest, ...
    %                     reconstructionDir, false, referenceImage, false, false);
    % 
    %                 %[ Store results
    %                 [~, registered.ms, registered.test] = load_niftyreg_results(...
    %                     bmsRegManifest, false, true, true) ;
    % 
    %                 registered.ms = concatenateCellTimeseries2Matrix(registered.ms);
    %                 registered.test = concatenateCellTimeseries2Matrix(registered.test);
    % 
    %                 plotdir = fullfile(fileparts(bmsRegManifest.registeredMotionStateFiles{1}),'plots');
    %                 if ~exist(plotdir,'dir'), mkdir(plotdir); end
    %                 % outstats = analyze_registered_MS_volumes( registered.ms, plotdir );
    % 
    %                 outgifname = fullfile(plotdir,'registered-MS-movie.gif');
    %                 createTimeseriesSlicesGif(registered.ms, 136, 94, 48, outgifname);
    % 
    %             end
    %         end
    %     end
    % end

    % if doPostBreathingMoco
    %     % Breathing motion correction 
    %     mcProjManifest = test_radial_vibe_28_mc_proj_imFIAT_with_pad(scanDir, projManifest, ...
    %         motionSignalManifest.zSI, currentTransforms, iReferenceMotionState, ...
    %         reconstructionDir, referenceImage, templateFile, smin_breath, smax_breath, ...
    %         stride_breath, false); % MC
    % end

    
end



    
