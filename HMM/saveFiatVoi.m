function saveFiatVoi(A, voiName, voiColor, fileName)

fid = fopen(fileName,'w');

[dstDir, dstFileNameRoot] = fileparts(fileName);

description = swapbytes(uint16(sprintf('%s:%s', voiName, fullfile(dstDir, dstFileNameRoot))));

fileHeader = swapbytes(uint16([0    513      0    2*numel(description)]));

fwrite(fid, fileHeader, 'uint16');

fwrite(fid, description, 'uint16');

fwrite(fid, uint8(255), 'uint8');

fwrite(fid, uint8(voiColor), 'uint8');

arraySize = cell(3, 1);
[arraySize{:}] = size(A);
arraySize = cell2mat(arraySize);
arraySize = typecast(swapbytes(uint16(arraySize)),'uint8');
fwrite(fid, arraySize, 'uint8');

iiSlices = fliplr(reshape(find(any(any(A,1),2)),1,[]));

for iSlice = iiSlices
    writeVoiSlice(fid, A, iSlice);
end

fwrite(fid, swapbytes(uint16(size(A,3))), 'uint16');

fclose(fid);

end

function writeVoiSlice(fid, A, iSlice)

fwrite(fid, swapbytes(uint16(iSlice - 1)), 'uint16');

iiMatlabColumns = find(any(A(:,:,iSlice),1));
iiMatlabRows = find(any(A(:,:,iSlice),2));

iStartRow = min(iiMatlabRows) - 1;
iEndRow = max(iiMatlabRows) - 1;
iStartColumn = min(iiMatlabColumns) - 1;
iEndColumn = max(iiMatlabColumns) - 1;

header = swapbytes(uint16([iStartRow iEndRow iStartColumn iEndColumn 0]));

fwrite(fid, header, 'uint16');

A1 = (A(:,:,iSlice) ~= 0);
A2 = imerode(A1,[0 1 0; 1 1 1; 0 1 0]);

fwrite(fid, uint8(A1 + A2), 'uint8');

end
