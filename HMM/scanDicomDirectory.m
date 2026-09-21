function t = scanDicomDirectory(srcDir, tagNames)

%files = dir(srcDir);

files = dirdicom(srcDir);

files = files(~strcmp({files(:).name},{'..'})); % isdir does not work for .. if path is prefixed by \\?\
files = files(~strcmp({files(:).name},{'.'})); % isdir does not work for . if path is prefixed by \\?\
files = files(~[files.isdir]);

nFiles = numel(files);


%tagNames = {'SeriesDescription', 'StudyInstanceUID', 'SeriesInstanceUID', 'PatientID', 'SOPInstanceUID'};
%tagCodes = {'(0008,103E)','(0020,000D)','(0020,000E)','(0010,0020)'};
%tagCodes = cellfun(@dicomHexPairToNumber, tagCodes);
[~, ~, DLN] = extractDicomDictionaryFromStandardNumeric();
tagCodes = cellfun(@(x)DLN(x), tagNames);
%tagCodes2 = arrayfun(@(x)dicomNumberToHexPair(x), tagCodes, 'UniformOutput', false); % Check

nTags = numel(tagCodes);

t = cell2table(cell(nFiles, numel(tagNames) + 1), 'VariableNames', [tagNames, {'fileName'}]);

m = false(nFiles,1);

% pg = Progressor('Scanning DICOM files...');

parfor iFile = 1:nFiles
%     try
        fileName = fullfile(srcDir,files(iFile).name);
        
        parsedTags = test_extract_dicom_tag(fileName, tagCodes);
        
        for iTag = 1:nTags
            switch tagNames{iTag}
                case 'PatientID'
                    t(iFile,:).(tagNames{iTag}){1} = normalizePatientId(char(parsedTags.body(([parsedTags.body.t] == tagCodes(iTag))).v'));
                otherwise
                    jTag = find([parsedTags.body.t] == tagCodes(iTag));
                    if numel(jTag) ~= 1
                        t(iFile,:).(tagNames{iTag}){1} = '';
                    else
                        t(iFile,:).(tagNames{iTag}){1} = char(parsedTags.body(([parsedTags.body.t] == tagCodes(iTag))).v');
                    end
            end
            %t{iFile, iTag} = char(parsedTags.body(([parsedTags.body.t] == tagCodes(iTag))).v');
        end  
        
        %dicomHexPairToNumber('(0002,0003)')
%         t.MediaStorageSOPInstanceUID{iFile} = strtrim(char(parsedTags.header([parsedTags.header.t] == 196610).v'));
        t(iFile,:).fileName{1} = fileName;
        m(iFile) = true;
%     catch
%         keyboard
%     end
    
%     if (mod(iFile, 15) == 1)
%         pg.setProgress(iFile/nFiles); 
%         pg.setDescription(sprintf('Scanning DICOM files... (%.0f/%.0f)',iFile,nFiles));
%     end
end

t = t(m,:);