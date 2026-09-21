function allDicomFiles = scanDicomDirectoryRecursively(srcDir, tagNames)

%fprintf('%s\n', srcDir);

dirs = dir(srcDir);
dirs = dirs(~strcmp({dirs(:).name},{'..'}));
dirs = dirs(~strcmp({dirs(:).name},{'.'}));
dirs = dirs([dirs.isdir]);

%dicomFiles = dirdicom(srcDir);
dicomFiles = scanDicomDirectory(srcDir, tagNames);
%t = scanDicomDirectory(srcDir, tagNames);

% for iFile = 1:numel(dicomFiles)
% %     fprintf('%s\n',fullfile(srcDir, dicomFiles(iFile).name));
%     %data = dicominfo(fullfile(srcDir, dicomFiles(iFile).name));
%     %fprintf('%s\n',data.SOPInstanceUID);
% end

allDicomFiles = cell(1,numel(dirs) + 1);
allDicomFiles{1} = dicomFiles;

parfor iDir = 1:numel(dirs)
    allDicomFiles{1 + iDir} = scanDicomDirectoryRecursively(fullfile(srcDir, dirs(iDir).name), tagNames);
end

allDicomFiles = cat(1, allDicomFiles{:});
