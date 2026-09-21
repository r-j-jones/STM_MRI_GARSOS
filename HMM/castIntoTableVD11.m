function mdh = castIntoTableVD11(mdhArray, mdhFields)

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
        mdhFields(iMdhField).data = typecast(mdhFields(iMdhField).data, mdhFields(iMdhField).type);
        nBlock = numel(mdhFields(iMdhField).data)/nMdh;
        mdhFields(iMdhField).data = reshape(mdhFields(iMdhField).data,[nBlock nMdh]);
        mdhFields(iMdhField).data = mdhFields(iMdhField).data.';
    else
        mdhFields(iMdhField).data = castIntoTableVD11(mdhFields(iMdhField).data, mdhFields(iMdhField).type);
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