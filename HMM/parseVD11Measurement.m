function m = parseVD11Measurement(measurementData, pStart, pEnd, varargin)

headerOnly = false;
calibrationOnly = false;

for i = 1:2:length(varargin)
    switch(varargin{i})
        case 'HeaderOnly'
            headerOnly = varargin{i+1};
        case 'CalibrationOnly'
            calibrationOnly = varargin{i+1};
        otherwise
            error('Unrecognized option: %s',varargin{i})
    end
end

headerLength = double(typecast(measurementData(pStart + (0:3)),'uint32'));

nBuffers = double(typecast(measurementData(pStart + (4:7)),'uint32'));

p = pStart + double(uint32(9)) - 1;

% pg = Progressor('Reading buffers.');
fprintf('Parsing buffers ('); tic;
tic;

for iBuffer = 1:nBuffers
    bufferName = measurementData(p + (0:31));
    bufferNameLength = (find(bufferName == 0, 1, 'first') - 1);
    bufferName = char(bufferName(1:bufferNameLength)');
    %pg.setDescription(sprintf('Reading buffers (%s).', bufferName));
    fprintf('%s', bufferName);
    p = p + bufferNameLength + 1;
    bufferLength = double(typecast(measurementData(p + (0:3)),'uint32'));
    p = p + 4;
    m.buffer.(bufferName) = char(measurementData(p + (0:bufferLength-2))');
    m.parsed.(bufferName) = parseXProtocol(m.buffer.(bufferName));
    m.param.(bufferName) = cell(size(m.parsed.(bufferName)));
    for ix = 1:numel(m.parsed.(bufferName))
        m.param.(bufferName){ix} = extractBufferParameters(m.parsed.(bufferName){ix});
    end
    p = p + bufferLength;
    %    pg.setProgress(double(iBuffer)/double(nBuffers));
    if iBuffer ~= nBuffers
        fprintf(', ');
    else
        fprintf(').');
        fprintf(' (%.2fs)\n',toc);
    end
end

m.header = m.buffer;

m.alignment = measurementData(p:pStart + headerLength - 1);

measurementIsCalibration = strcmp(m.param.Meas{1}.MEAS.tProtocolName,'AdjCoilSens');

if ~headerOnly && (~calibrationOnly || measurementIsCalibration)
    
    p = pStart + headerLength;
    
    nMdh = 0;
    
    fileSize = pEnd - pStart + 1;
    dataSize = fileSize - headerLength;
    
    % Read samples
    
    maxNMdh = floor(dataSize/192);
    
    pMdh = zeros(maxNMdh, 1, 'double');
    rawDataLength = zeros(maxNMdh, 1, 'double');
    rawDataNCoils = zeros(maxNMdh, 1, 'double');
    
    endOfScan = false;
    
    fprintf('Scanning data.'); tic;
    
    while ~endOfScan
        
        q0 = measurementData(p + 40 + (0:7));
        evalInfoMask = typecast(q0, 'uint64');
        
        isSyncData = logical(bitshift(bitand(evalInfoMask,bitshift(uint64(1), 5)),-5));
        
        if isSyncData
            ulDMALength = double(typecast(measurementData(p + (0:3)),'uint32'));
            p = p + ulDMALength;
        else
            
            nMdh = nMdh + 1;
            pMdh(nMdh) = p;
            
            q1 = measurementData(p + 48 + (0:1));
            q2 = measurementData(p + 48 + (2:3));
            rawDataLength(nMdh) = double(typecast(q1, 'uint16'));
            rawDataNCoils(nMdh) = double(typecast(q2, 'uint16'));
            %     assert(rawDataLength == 384,'Oh nooo!');
            endOfScan = logical(bitand(measurementData(p + 40), 1));
            p = p + 192 + (8*rawDataLength(nMdh) + 32)*rawDataNCoils(nMdh);
            %     if mod(nMdh, 100) == 1
            %         pg.setProgress(double(p - headerLength)/double(dataSize));
            %     end
            
        end
        
    end
    
    pMdh = pMdh(1:nMdh);
    rawDataLength = rawDataLength(1:nMdh);
    rawDataNCoils = rawDataNCoils(1:nMdh);
    
    mdhArray = zeros(192, nMdh, 'uint8');
    % coilHeaderArray = cell(1, nMdh);
    channelHeaderArray = cell(1, nMdh);
    dataArray = cell(1, nMdh);
    
    fprintf(' (%.2fs)\n',toc);
    fprintf('Segmenting data.'); tic;
    
    for iMdh = 1:nMdh
        mdhArray(:,iMdh) = measurementData(pMdh(iMdh) + (0:191));
        lineDataPartitionSizes = [32, 8*rawDataLength(iMdh)];
        lineBlockSize = [sum(lineDataPartitionSizes), rawDataNCoils(iMdh)];
        pSamplesBegin = pMdh(iMdh) + 192;
        pSamplesEnd = pSamplesBegin + prod(lineBlockSize) - 1;
        
        q1 = measurementData(pSamplesBegin:pSamplesEnd);
        
        q2 = reshape(q1, lineBlockSize);
        
        q3 = mat2cell(q2, lineDataPartitionSizes, rawDataNCoils(iMdh));
        
        %     q = measurementData(pSamplesBegin:pSamplesEnd);
        %     q = typecast(q, 'single');
        % %     q1 = dataArray{iMdh}(1:2:end);
        % %     q2 = dataArray{iMdh}(2:2:end);
        %     q = reshape(q,2,[]);
        %     dataArray{iMdh} = complex(q(1,:), q(2,:));
        % %     if mod(iMdh, 1000) == 1
        % %         pg.setProgress(iMdh/nMdh);
        % %     end
        channelHeaderArray{iMdh} = q3{1};
        dataArray{iMdh} = q3{2};
    end
    
    fprintf(' (%.2fs)\n',toc);
    fprintf('Casting data.'); tic;
    
    %m.channelHeaderArray = channelHeaderArray;
    %m.dataArray = dataArray;
    
    mdhFields = getMdhStructureVD11();
    
    mdhTable = castIntoTableVD11(mdhArray, mdhFields);
    
    channelHeaderFields = getChannelHeaderStructureVD11();
    
    channelHeaderArray = cell2mat(channelHeaderArray);
    channelHeaderArray = castIntoTableVD11(channelHeaderArray, channelHeaderFields);
    
    dataArray = cell2mat(cellfun(@(x)x(:), dataArray, 'UniformOutput', false)');
    dataArray = typecast(dataArray,'single');
    dataArray = complex(dataArray(1:2:end), dataArray(2:2:end));
    dataArrayLengths = cell2mat(arrayfun(@(x,y)repmat(x, [y 1]), double(mdhTable.ushSamplesInScan), double(mdhTable.ushUsedChannels), 'UniformOutput', false));
    dataArray = mat2cell(dataArray, dataArrayLengths);
    
    channelHeaderArray = [channelHeaderArray, table(dataArray, 'VariableNames', {'rawdata'})];
    
    channelHeaderArray = mat2cell(channelHeaderArray, double(mdhTable.ushUsedChannels), size(channelHeaderArray, 2));
    
    mdhTable = [mdhTable, table(channelHeaderArray,'VariableNames',{'channelHeaders'})];
    
    m.mdh = mdhTable;
    
    % Cast samples
    
    %m = rmfield(m,{'parsed','alignment','buffer'});
    %m = rmfield(m,{'alignment','buffer','dataArray'});
    
    % m.pMdh = pMdh;
    % m.rawDataLength = rawDataLength;
    % m.rawDataNCoils = rawDataNCoils;
    % m.dataLength = numel(measurementData);
    
    fprintf(' (%.2fs)\n',toc);
    
end