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
        metaManifest, outputDir, referenceImage, opts.forceUpdate);


    disp('done getting stm data...');


    %[ Save the results
    [garDataDir,~,~] = fileparts(fmanifests.proj);
    outputMatFile = fullfile(garDataDir,'test_STM_GARSOS_data-prep_09-20-2026.mat');
    load(outputMatFile,'inputOptions');


    


    %% STM computation and analysis

    imageSize = projManifest.sequences.timeSeries.imageSize;
    N1 = imageSize(1);
    N2 = imageSize(2);

    kCalVolume = projManifest.sequences.kCal;
    [N1_cal, N2_cal, Nz, Nt] = size(kCalVolume);  % N1 x N2 : image dimensions
                                    % Nc      : number of coils
                                    % Nt      : number of time frames


    %% STM computation parameters

    dim_sens = [N1, N2];                  % Desired dimensions for the computed spatiotemporal maps.
    
    L = 7;                                % Number of temporal basis functions to be computed.
    
    tau      = 3;                         % Kernel radius. Default: 3
    
    threshold = 0.035;                     % Threshold for C-matrix singular values. Default: 0.05
                                          % Note: In this example we don't use the default value.
    
    M = 20;                               % Number of iterations for Orthogonal Iteration. Default: 30
                                          % Note: In this example we use a smaller value
                                          % to speed up the calculations.
    
    interp_zp = 192;                       % Amount of zero-padding to create the low-resolution grid 
                                          % if FFT-interpolation is used. Default: 24
    
    gauss_win_param = 100;                % Parameter for the Gaussian apodizing window used to 
                                          % generate the low-resolution image in the FFT-based 
                                          % interpolation approach. This is the reciprocal of the 
                                          % standard deviation of the Gaussian window. Default: 100
    
    sketch_dim = 300;                     % Dimension of the sketch matrix used to calculate a
                                          % basis for the nullspace of the C matrix using a sketched SVD. 
                                          % Default: 500. Note: In this example we use a smaller value
                                          % to speed up the calculations.
    
    visualize_C_matrix_sv = 1;            % Binary variable. 1 = Singular values of the C matrix are displayed.
                                          % Default: 0. 
                                          % Note: In this example we set it to 1 to visualize the singular values
                                          % of the C matrix. If sketched_SVD = 1 and if the curve of the singular values flattens out,
                                          % it suggests that the sketch dimension is appropriate for the data.
                                          
    %% PISCO techniques
    
    % The following techniques are used if the corresponding binary variable is equal to 1
    
    kernel_shape = 1;                     % Binary variable. 1 = ellipsoidal shape is adopted for 
                                          % the calculation of kernels (instead of rectangular shape).
                                          % Default: 1
    
    FFT_nullspace_C_calculation = 1;      % Binary variable. 1 = FFT-based calculation of nullspace 
                                          % vectors of C by calculating C'*C directly (instead of 
                                          % calculating C first). Default: 1
    
    sketched_SVD = 1;                     % Binary variable. 1 = sketched SVD is used to calculate 
                                          % a basis for the nullspace of the C matrix (instead of 
                                          % calculating the nullspace vectors directly and then the 
                                          % basis). Default: 1
    
    OrthogonalIteration_G_nullspace_vectors = 0; % Binary variable. 1 = Orthogonal Iteration approach is 
                                            % used to find nullspace vectors of the G matrices 
                                            % (instead of using SVD). Default: 1
    
    FFT_interpolation = 1;                % Binary variable. 1 = sensitivity maps are calculated on 
                                          % a small spatial grid and then interpolated to a grid with 
                                          % nominal dimensions using an FFT-approach. Default: 1
    
    verbose = 1;                          % Binary variable. 1 = PISCO information is displayed. 
                                          % Default: 1
    
    %% PISCO estimation
    
    if isempty(which('STM_computation'))
        error(['The function STM_computation.m is not found in your MATLAB path. ' ...
               'Please ensure that all required files are available and added to the path.']);
    end

    [garDataDir,~,~] = fileparts(fmanifests.proj);
    plotDir = fullfile(garDataDir,'stm-plots');
    if ~exist(plotDir,'dir')
        mkdir(plotDir);
    end


    
    kCal = squeeze(kCalVolume(:,:,50,:));
    
    t_stm = tic;
    
    [ST_maps, eigenValues] = STM_computation( ...
        kCal, ...
        dim_sens, ...                          % Data and output size
        L, ...
        'tau', tau, ...
        'threshold', threshold, ...
        'kernel_shape', kernel_shape, ...            % Kernel and threshold parameters
        'FFT_nullspace_C_calculation', FFT_nullspace_C_calculation, ...             % FFT nullspace calculation flag
        'OrthogonalIteration_G_nullspace_vectors', OrthogonalIteration_G_nullspace_vectors, ...      % Orthogonal Iteration flag
        'M', M, ...
        'FFT_interpolation', FFT_interpolation, ...
        'interp_zp', interp_zp, ...
        'gauss_win_param', gauss_win_param, ... % Interpolation params
        'sketched_SVD', sketched_SVD, ...
        'sketch_dim', sketch_dim, ...
        'visualize_C_matrix_sv', visualize_C_matrix_sv, ... % SVD/sketching params
        'verbose', verbose ...                                  % Verbosity
    );
    
    disp(['Time for STM computation: ' num2str(toc(t_stm)) ' seconds']);
    disp('=======================');

    if visualize_C_matrix_sv
        f=gcf;
        outPlotPath = fullfile(plotDir,'C_matrix_SV.png');
        print(f,outPlotPath,'-dpng','-r300');
    end
    
    f2 = figure; 
    imagesc(utils.mdisp(abs(eigenValues)));
    axis tight;
    axis image;
    colorbar; 
    colormap gray;
    clim([0 1])
    title('Eigenvalues of G matrices (normalized)');
    f2.Position = [570 82 989 407];
    outPlotPath = fullfile(plotDir,'G_matrix_eigenvals.png');
    print(f2,outPlotPath,'-dpng','-r300');


    %% Save ST_maps and eigenValues

    outputMatFile = fullfile(garDataDir,'reconed-ST_maps.mat');
    save(outputMatFile,'ST_maps','eigenValues');


    % %% Check example STM for dims
    % 
    % cd('/RadOnc-MRI1/Student_Folder/rjones/STM/STM_MRI_GARSOS');
    % example_STM_recon_2D_multichannel;
    % 
    % 
    % % Orig ST_maps = 128x128x40x4
    % % Orig kCal = 128     8    40
    % % Orig interp_zp = 100
    % 
    %% STM reconstruction
    
    N = N1*N2;
    
    % Operators to apply spatiotemporal maps
    
    Sv = @(x) utils.vect(sum(ST_maps.*repmat(reshape(x, [N1 N2 1 L]), [1 1 Nt 1]), 4)); % Operators to apply temporal maps
    Sv_h = @(x) utils.vect(sum(conj(ST_maps).*repmat(reshape(x, [N1 N2 Nt 1]), [1 1 1 L]), 3));
    Sv_h_Sv = @(x) Sv_h(Sv(x));
    
    % Operators to apply sensitivity maps
    
    Fv = @(x) utils.vect(repmat(reshape(x, [N1 N2 1 Nt]), [1 1 Nc 1]) .* repmat(sense_maps, [1 1 1 Nt]));
    Fv_h = @(x) utils.vect(sum(conj(sense_maps).*reshape(x, [N1 N2 Nc Nt]), 3));
    
    % Forward system operator
    
    A = @(x) utils.vect(kmask.*utils.ft2(reshape(Fv(Sv(x)), [N1 N2 Nc Nt])));
    
    % Adjoint system operator
    
    Ah = @(x) Sv_h(Fv_h(utils.vect(utils.ift2(kmask.*reshape(x, [N1 N2 Nc Nt])))));
    
    % Composition operator
    
    AhA =  @(x) Ah(A(x));
    
    %% CG-STM reconstruction
    
    % This is a simple reconstruction with no regularization
    
    disp('Starting CG-STM reconstruction...');
    disp('=======================');
    
    t_stm_recon = tic;
    
    [z, ~] = pcg(AhA, Ah(kdata_under), 1e-6, 50); % Reconstructed STM spatial coefficients
    
    sense_recon_stm = reshape(Sv(z), [N1 N2 Nt]);
    
    disp(['Time for CG-STM reconstruction: ' num2str(toc(t_stm_recon)) ' seconds']);
    disp('=======================');
    
    % NRMSE CG-STM reconstruction
    
    NRMSE_stm = norm(idata_gt_sc(:) - sense_recon_stm(:))/norm(idata_gt_sc(:));
    
    disp(['NRMSE CG-STM reconstruction: ' num2str(NRMSE_stm)]);
    disp('=======================');
    
    % Visualization of CG-STM reconstruction
    
    figure;
    imagesc(utils.mdisp(abs(sense_recon_stm))); 
    colormap gray;
    axis tight;
    axis image;
    axis off;
    title('CG-STM reconstruction - all frames');
    
    %% CG-STM + Tikhonov reconstruction
    
    % Reconstruction with a simple Tikhonov regularization
    
    lambda_tik = 0.001; % Tikhonov regularization parameter
    
    AtikA = @(x) AhA(x) + lambda_tik*x;
    
    disp('Starting CG-STM + Tikhonov reconstruction...');
    disp('=======================');
    
    t_stm_recon_tik = tic;
    
    [z_tik, ~] = pcg(AtikA, Ah(kdata_under), 1e-6, 50); % Reconstructed STM spatial coefficients
    
    sense_recon_stm_tik = reshape(Sv(z_tik), [N1 N2 Nt]);
    
    disp(['Time for CG-STM + Tikhonov reconstruction: ' num2str(toc(t_stm_recon_tik)) ' seconds']);
    disp('=======================');
    
    % NRMSE CG-STM + Tikhonov reconstruction
    
    NRMSE_stm_tik = norm(idata_gt_sc(:) - sense_recon_stm_tik(:))/norm(idata_gt_sc(:));
    
    disp(['NRMSE CG-STM + Tikhonov reconstruction: ' num2str(NRMSE_stm_tik)]);
    disp('=======================');
    
    % Visualization of CG-STM + Tikhonov reconstruction
    
    figure;
    imagesc(utils.mdisp(abs(sense_recon_stm_tik))); 
    colormap gray;
    axis tight;
    axis image;
    axis off;
    title('CG-STM + Tikhonov reconstruction - all frames');
    





    if 0 == 1

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

    end

    
end



    
