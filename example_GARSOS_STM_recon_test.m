%% Script to test GARSOS STM

inputPath = '/mnt/ibrixfs04-Kspace/motion_patients/P0021/E02';
outputPath = '/RadOnc-MRI1/Student_Folder/rjones/STM/InitialTest_09-20-2026';
opts = [];
opts.nbSpokesPerFrame = 15;

% Load reference (vendor-reconed) dicom volume
referenceImageDir = fullfile(outputPath,'reference-dcm-vol');
if ~exist(referenceImageDir,'dir'), mkdir(referenceImageDir); end
referenceImageFile = fullfile(referenceImageDir,'referenceVolumeA.mat');
load(referenceImageFile,'referenceVolumeA');
opts.referenceVolumeA = referenceVolumeA;

cd('/RadOnc-MRI1/Student_Folder/rjones/STM/STM_MRI_GARSOS/HMM');

result = main_STM_GARSOS(inputPath, outputPath, opts);

