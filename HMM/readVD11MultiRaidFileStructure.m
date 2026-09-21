function s = readVD11MultiRaidFileStructure(raidFile, varargin)

%raidFile = 'D:\adamroot\data\MRraw3\VD11test\meas.dat';

fprintf('Mapping file into memory.'); tic;

mapFile = memmapfile(raidFile);

fileData = mapFile.Data;

s.hdSize_ = typecast(fileData(1:4),'uint32');
s.count_ = typecast(fileData(5:8),'uint32');

mrParcRaidFileEntryFields = getMrParcRaidFileEntryStructure();

s.mrParcRaidFileEntries = castIntoTable(reshape(fileData(8 + (1:152*64)),[152 64]),mrParcRaidFileEntryFields);

s.m = cell(1, s.count_);

fprintf(' (%.2fs)\n',toc);

for iMeasurement = 1:s.count_
    fprintf('Reading measurement #%.0f of %.0f.\n', iMeasurement, s.count_);
    pStart = s.mrParcRaidFileEntries(iMeasurement,:).off_ + 1;
    pEnd = pStart + s.mrParcRaidFileEntries(iMeasurement,:).len_ - 1;
    s.m{iMeasurement} = parseVD11Measurement(fileData, double(pStart), double(pEnd),varargin{:});
end

fprintf('Releasing file.'); tic;

clear fileData;

clear mapFile;

fprintf(' (%.2fs)\n',toc);

% function m = parseVD11Measurement(measurementData) %#ok
% 
% m = [];

function mdh = castIntoTable(mdhArray, mdhFields)

blockLength = sum([mdhFields.length]);

% fprintf('blockLength = %.0f\n', blockLength);

nMdh = numel(mdhArray)/sum([mdhFields.length]);

mdhArray = reshape(mdhArray, [blockLength, nMdh]);

mdhArray = mat2cell(mdhArray,[mdhFields.length],nMdh);

% pg = Progressor('Cast into fields.');

for iMdhField = 1:numel(mdhFields)
    mdhFields(iMdhField).data = mdhArray{iMdhField};
    if ischar(mdhFields(iMdhField).type)
        mdhFields(iMdhField).data = reshape(mdhFields(iMdhField).data,1,mdhFields(iMdhField).length*nMdh);
        if strcmp(mdhFields(iMdhField).type, 'char')
            mdhFields(iMdhField).data = char(typecast(mdhFields(iMdhField).data, 'uint8'));
        else
            mdhFields(iMdhField).data = typecast(mdhFields(iMdhField).data, mdhFields(iMdhField).type);
        end
        nBlock = numel(mdhFields(iMdhField).data)/nMdh;
        mdhFields(iMdhField).data = reshape(mdhFields(iMdhField).data,[nBlock nMdh]);
        mdhFields(iMdhField).data = mdhFields(iMdhField).data.';
    else
        mdhFields(iMdhField).data = castIntoTable(mdhFields(iMdhField).data, mdhFields(iMdhField).type);
    end
    %%nBlock = numel(mdhFields(iMdhField).data)/nMdh;
    %%mdhFields(iMdhField).data = reshape(mdhFields(iMdhField).data,[nBlock nMdh]);
    %     mdhFields(iMdhField).data = mat2cell(mdhFields(iMdhField).data,nBlock,ones(1,nMdh));
    %mdhFields(iMdhField).data = num2cell(mdhFields(iMdhField).data,1);
    %%mdhFields(iMdhField).data = mdhFields(iMdhField).data.';
    
%     pg.setProgress(iMdhField/numel(mdhFields));
end
    
mdhStructConstructor = cell(2,numel(mdhFields));

mdhStructConstructor(1,:) = {mdhFields.name};
mdhStructConstructor(2,:) = {mdhFields.data};

%mdh = struct(mdhStructConstructor{2,:},'VariableNames',mdhStructConstructor(1,:));
mdh = table(mdhStructConstructor{2,:},'VariableNames',mdhStructConstructor(1,:));

function mrParcRaidFileEntryFields = getMrParcRaidFileEntryStructure()

mrParcRaidFileEntryFields(1).name = 'measId_';
mrParcRaidFileEntryFields(1).type = 'uint32';
mrParcRaidFileEntryFields(1).length = 4;

mrParcRaidFileEntryFields(2).name = 'fileId_';
mrParcRaidFileEntryFields(2).type = 'uint32';
mrParcRaidFileEntryFields(2).length = 4;

mrParcRaidFileEntryFields(3).name = 'off_';
mrParcRaidFileEntryFields(3).type = 'uint64';
mrParcRaidFileEntryFields(3).length = 8;

mrParcRaidFileEntryFields(4).name = 'len_';
mrParcRaidFileEntryFields(4).type = 'uint64';
mrParcRaidFileEntryFields(4).length = 8;

mrParcRaidFileEntryFields(5).name = 'patName_';
mrParcRaidFileEntryFields(5).type = 'char';
mrParcRaidFileEntryFields(5).length = 64;

mrParcRaidFileEntryFields(6).name = 'protName_';
mrParcRaidFileEntryFields(6).type = 'char';
mrParcRaidFileEntryFields(6).length = 64;
