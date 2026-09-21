function [arrayVolume, outputGridPoints] = loadArrayVolumeFromDicomSeriesInstanceUid3(dicomFiles, indexNames)

setAugmentedDicomDictionary();

%metaData = num2cell(cellfun(@(x)dicominfo(x),dicomFiles));
%metaData = cellfun(@(x)dicominfo(x), dicomFiles, 'UniformOutput', false);

nDicomFiles = numel(dicomFiles);

metaData = cell(1, nDicomFiles);

if ~is_in_parallel()
    
    p = gcp;
    
    pctRunOnAll('setAugmentedDicomDictionary()');
    
end

parfor iDicomFile = 1:nDicomFiles
    metaData{iDicomFile} = dicominfo(dicomFiles{iDicomFile});
end

% TODO:
indexTable = getIndexVariablesFromDicomHeaders(metaData, indexNames);

% load the instance num, rows and columns attribute in the dicom file meta info
nInstances = numel(metaData);

nRows = metaData{1}.Rows;
nColumns = metaData{1}.Columns;

R = reshape(metaData{1}.ImageOrientationPatient,[3 2]);
R = [R, cross(R(:, 1), R(:, 2))];

r0 = cell2mat(cellfun(@(x)x.ImagePositionPatient, metaData, 'UniformOutput', false));

d = R\r0;

indexNames = [{'SLICE_LOCATION'}, indexNames];

indexTable = cat(2, table(d(3,:)','VariableNames',indexNames(1)), indexTable);

%[~, iiSlice] = sort(d(3, :), 'ascend');

[indexTable, iiInstance] = sortrows(indexTable, fliplr(indexNames),'ascend');

metaData = metaData(1, iiInstance);

d = d(:, iiInstance);
r0 = r0(:, iiInstance);

gridPoints = cellfun(@(x)sort(unique(indexTable.(x))), indexNames, 'UniformOutput', false);

outputGridPoints = gridPoints(2:end);

arraySize = [nColumns, nRows, cellfun(@(x)numel(x), gridPoints)];

assert(prod(arraySize(3:end)) == nInstances, 'Cannot reshape %s.', dicomFiles{1});

A = zeros(arraySize);

%pg = Progressor('Loading voxel data.');

parfor iInstance = 1:nInstances
    A(:,:,iInstance) = transpose(dicomread(metaData{iInstance}));
    %pg.setProgress(iInstance/nInstances);
end

sliceSpacing = diff(gridPoints{1});

nUniqueSliceSpacings = numel(unique(round(sliceSpacing/0.001)));

assert(nUniqueSliceSpacings <= 1, 'Inconsistent slice spacing.');

if nUniqueSliceSpacings == 1
    sliceSpacing = sliceSpacing(1);
elseif nUniqueSliceSpacings == 0
    if isfield(metaData{1}, 'SliceThickness')
        sliceSpacing = metaData{1}.SliceThickness;
    else
        sliceSpacing = 1;
    end
end


voxelSpacing = [metaData{1}.PixelSpacing(2); metaData{1}.PixelSpacing(1); sliceSpacing];

arrayVolume = ArrayVolume(A, R, voxelSpacing, r0(:,1), reshape(metaData, [arraySize(3:end) 1 1]));






