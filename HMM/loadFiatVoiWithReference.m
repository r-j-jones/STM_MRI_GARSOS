function [roiVolume, voiName, voiColor] = loadFiatVoiWithReference(voiFile, refVolume)

[A, voiName, voiColor] = loadFiatVoi(voiFile);

refSize = cell(1, 4);
[refSize{:}] = size(refVolume.A);
refSize = cell2mat(refSize);

roiSize = cell(1, 4);
[roiSize{:}] = size(A);
roiSize = cell2mat(roiSize);

assert(all(refSize(1:3) == roiSize(1:3)));

roiVolume = ArrayVolume(A, refVolume.R, refVolume.v, refVolume.r0, []);
