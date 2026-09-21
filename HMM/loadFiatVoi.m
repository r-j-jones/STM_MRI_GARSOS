function [A, voiName, voiColor, originalFileName, B] = loadFiatVoi(srcFile)

assert(exist(srcFile,'file') == 2, 'The file does not exist. (%s)', srcFile);

fid = fopen(srcFile);

t = fread(fid,'uint8=>uint8')';

fclose(fid);

fileHeader = swapbytes(typecast(t(1:8),'uint16'));

stringLength = double(fileHeader(4));

fileNameEndMark = 8 + stringLength + 1;

description = char(swapbytes(typecast(t(9:fileNameEndMark - 1),'uint16')));

iSeparator = find(description == ':');

voiName = description(1:iSeparator - 1);
originalFileName = description(iSeparator + 1:end);

voiColor = double(t(fileNameEndMark + (1:3)));

arraySize = double(swapbytes(typecast(t(fileNameEndMark + 3 + (1:6)),'uint16')));
%Order: DICOM iDicomColumn,iDicomRow,iDicomSlice <=> MATLAB iMatlabRow,iMatlabCol,iMatlabSlice

A = zeros(arraySize);
B = zeros(arraySize);

p = fileNameEndMark + 3 + 6 + 1;

while p < numel(t)
    [pnew, iSlice, sliceData] = readVoiSlice(t, p, arraySize);
    p = pnew;
    if iSlice <= arraySize(3)
        A(:,:,iSlice) = (mod(sliceData, 3) ~= 0);
        B(:,:,iSlice) = sliceData;
    end
end

end

function [p, iSlice, sliceData] = readVoiSlice(t, p0, arraySize)

iSlice = double(swapbytes(typecast(t(p0 + (0:1)), 'uint16'))) + 1;

if iSlice <= arraySize(3)
    
    header = swapbytes(typecast(t(p0 + (2:11)), 'uint16'));
    
    p = p0 + 12;
    
    sliceData = reshape(t(p + (0:prod(arraySize(1:2)) - 1)), arraySize(1:2));
    
    p = p + prod(arraySize(1:2));
    
else
    p = p0 + 2;
    sliceData = [];
end

end
