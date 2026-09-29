% Script to test GARSOS STM

inputPath = '/mnt/ibrixfs04-Kspace/motion_patients/P0021/E02';
outputPath = '/RadOnc-MRI1/Student_Folder/rjones/STM/InitialTest_09-20-2026';
opts = [];
opts.nbSpokesPerFrame = 15;

cd('/RadOnc-MRI1/Student_Folder/rjones/STM/STM_MRI_GARSOS/HMM');

result = main_STM_GARSOS(inputPath, outputPath, opts);

