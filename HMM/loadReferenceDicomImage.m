function [referenceImage, templateFile, otherOutput] = loadReferenceDicomImage( scanDir, parFlag )
% function [referenceImage, templateFile, otherOutput] = loadReferenceDicomImage( scanDir[, parFlag] )
%   parFlag, opt. (default = true)

if nargin<2 || isempty(parFlag), parFlag = true; end

fprintf('parFlag = %s\n',parFlag);

fprintf(' -loading reference dicom image..\n');

%%% Load reference Image here
if parFlag
    instanceTable = scanDicomDirectoryRecursively(fullfile(scanDir, 'dcm'), {'SeriesInstanceUID', 'StudyInstanceUID', 'PatientID', 'StudyDate', 'SeriesDescription', 'SeriesNumber', 'Modality'});
else
    disp('not using parfor loading...')
    instanceTable = scanDicomDirectoryRecursively_nopar_rj(fullfile(scanDir, 'dcm'), {'SeriesInstanceUID', 'StudyInstanceUID', 'PatientID', 'StudyDate', 'SeriesDescription', 'SeriesNumber', 'Modality'});
end
seriesTable = instanceTable(:, [6 1 5]);
seriesTable.SeriesNumber = cellfun(@(x)str2double(x), seriesTable.SeriesNumber);
seriesTable = unique(seriesTable,'rows');
subTable = instanceTable(strcmp(instanceTable.SeriesInstanceUID, seriesTable.SeriesInstanceUID{1}),:);
if parFlag
    referenceImage = loadArrayVolumeFromDicomSeriesInstanceUid3(subTable.fileName,{});
else
    disp('not using parfor loading...')
    referenceImage = loadArrayVolumeFromDicomSeriesInstanceUid3_nopar(subTable.fileName,{});
end
% save 'referenceImage.mat' referenceImage;
templateFile = string(instanceTable.fileName(1));

otherOutput.instanceTable = instanceTable;
otherOutput.seriesTable = seriesTable;
otherOutput.subTable = subTable;
