function metaData = getMetaDataForSeriesInstanceUid(seriesInstanceUid)

if ~iscell(seriesInstanceUid)
    assert(ischar(seriesInstanceUid));
    seriesInstanceUid = {seriesInstanceUid};
end

db = PostgreSqlDicomDatabase('dicomdb');

fileNames = cellfun(@(x)db.getConditionalDistinctTagValues('fileName','SeriesInstanceUID',x), seriesInstanceUid, 'UniformOutput', false);

fileNames = reshape(fileNames, [1 numel(fileNames)]);

fileNames = cellfun(@(x)reshape(x, 1, numel(x)), fileNames, 'UniformOutput', false);

fileNames = cat(2, fileNames{:});

metaData = cell(1, numel(fileNames));

nFiles = numel(fileNames);

if is_in_parallel()    
    setAugmentedDicomDictionary();
else
    p = gcp;
    pctRunOnAll('setAugmentedDicomDictionary()');
end

% setAugmentedDicomDictionary();

parfor iFile = 1:nFiles
    metaData{iFile} = dicominfo(fileNames{iFile});
end

% if iscell(seriesInstanceUid)
%     
%     %pg = Progressor('Loading meta data for series.');
%     
%     metaData = cell(size(seriesInstanceUid));
%     
%     nSeries = numel(seriesInstanceUid);
%     
%     if nSeries > 1
%     
%     parfor iSeries = 1:nSeries
%     %for iSeries = 1:nSeries
%         %pg.setDescription(sprintf('Loading meta data for series %.0f of %.0f.',iSeries, nSeries));
%         metaData{iSeries} = getMetaDataForSeriesInstanceUid(seriesInstanceUid{iSeries});
%         %pg.setProgress(iSeries/nSeries);
%     end
%     
%     metaData = cat(2, metaData{:});
%     
% else
%     
%     db = PostgreSqlDicomDatabase('dicomdb');
%     
%     fileNames = db.getConditionalDistinctTagValues('fileName','SeriesInstanceUID',seriesInstanceUid);
%     
%     metaData = cell(1, numel(fileNames));
%     
%     pg = Progressor('Reading DICOM meta data.');
%     
%     nFiles = numel(fileNames);
%     
%     for iFile = 1:nFiles;
%         metaData{iFile} = dicominfo(fileNames{iFile});
%         pg.setDescription(sprintf('Reading DICOM meta data (%.0f/%.0f).', iFile, nFiles));
%         pg.setProgress(iFile/nFiles);
%     end
%     
% end