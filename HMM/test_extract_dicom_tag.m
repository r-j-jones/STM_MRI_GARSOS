function data = test_extract_dicom_tag(srcFile, tagList)

memoryMap = memmapfile(srcFile);
m = memoryMap.data;

p = 1;

[p, ~] = readPreamble(p, m);
[p, x] = readSignature(p, m);

[D61, D71, DLN] = extractDicomDictionaryFromStandardNumeric();

recursionLevel = 0;

if strcmp(char(x'),'DICM')
    
    allowPrivate = false;
    [p, data.header] = parseDataElementsEarlyTermination(p, m, numel(m),D71, 'explicit', allowPrivate, recursionLevel, tagList);
    
    %transferSyntaxUid = char(data.header(([data.header.t] == dicomHexPairToNumber('(0002,0010)'))).v)';
    transferSyntaxUid = char(data.header(([data.header.t] == DLN('TransferSyntaxUID'))).v)';
    transferSyntaxUid(transferSyntaxUid == 0) = [];
    
    switch transferSyntaxUid
        case '1.2.840.10008.1.2'
            transferSyntax = 'implicit';
        case '1.2.840.10008.1.2.1'
            transferSyntax = 'explicit';
        otherwise
            error('Unknown transfer syntax uid: %s\n', transferSyntaxUid);
    end
    
else
    p = 1;
    groupNumber = typecast(m(1:2),'uint16');
    
    if groupNumber == uint16(2) || groupNumber == uint16(8)
        transferSyntax = 'implicit';
    else
        error('Could not read file  as DICOM.');
    end
end

allowPrivate = true;
[~, data.body] = parseDataElementsEarlyTermination(p, m, numel(m), D61, transferSyntax, allowPrivate, recursionLevel, tagList);


function [p, x] = readPreamble(p, m)

x = m(p + (0:127));
p = p + 128;

function [p, x] = readSignature(p, m)

x = m(p + (0:3));
%assert(all(char(x)' == 'DICM'),'DICOM prefix "DICM" not found.');
p = p + 4;

function [p, elementList] = parseDataElementsEarlyTermination(p, m, pEnd, dictionary, transferSyntax, allowPrivate, recursionLevel, tagList)

elementList = struct('t',{},'vr',{},'l',{},'v',{},'pStart',{},'pNextStart',[]);

x = 1;

tagsFound = false(size(tagList));

while ~all(tagsFound) && ~isempty(x) && p <= pEnd
    [p, x] = parseDataElement(p, m, pEnd, dictionary, transferSyntax, allowPrivate, recursionLevel, tagList);
    if ~isempty(x)
        elementList(end + 1) = x;
        tagsFound = tagsFound | (x.t == tagList);
    end
end

function [p, dataElement] = parseDataElement(p, m, pEnd, dictionary, transferSyntax, allowPrivate, recursionLevel, tagList)

pStart = p;

%padding = repmat(' ',[1 recursionLevel]);

% groupNumber = typecast(m(p + (0:1)),'uint16');
% elementNumber = typecast(m(p + (2:3)),'uint16');
% tag = double(typecast([groupNumber, elementNumber],'uint32'));

tag = double(typecast(m(p + (0:3)),'uint32'));

%tag = sprintf('(%0.4X,%0.4X)', groupNumber, elementNumber);

inDictionary = dictionary.isKey(tag);

%dicomHexPairToNumber('(FFFE,E00D)')

if (inDictionary || allowPrivate) && ~(tag == 3759013886) %~strcmp(tag, '(FFFE,E00D)')
    
    switch transferSyntax
        case 'implicit'
            if inDictionary
                c = dictionary(tag);
                vr = c{3};
            else
                vr = 'UN';
            end
            
            l = double(typecast(m(p + (4:7)),'uint32'));
            p = p + 8;
        case 'explicit'
            vr = char(m(p + (4:5))');
            
            if inDictionary
                c = dictionary(tag);
                if isempty(regexp(c{3}, vr, 'once'))
                    warning('Value representation (VR) does not match element.');
                end
                %assert(~isempty(regexp(c{3}, vr, 'once')),'Value representation (VR) does not match element.');
            end
            
            switch vr
                case {'OB', 'OW', 'OF', 'SQ', 'UR', 'UT', 'UN'}
                    l = double(typecast(m(p + (8:11)),'uint32'));
                    p = p + 12;
                otherwise
                    l = double(typecast(m(p + (6:7)),'uint16'));
                    p = p + 8;
            end
    end
    
%     fprintf('%s%s %8s %8.0f: ', padding, tag, vr, l);
    
%     if inDictionary
%         fprintf(' %s', c{1});
%     end
    
%     fprintf('\n');
    
    if l == 4294967295
        if strcmp(vr, 'UN')
            pSequenceStart = p;
            [p, ~] = parseSequenceItems(p, m, pEnd, dictionary, 'implicit', allowPrivate, recursionLevel + 1, tagList);
            v = m(pSequenceStart:p-1);
        else
            [p, v] = parseSequenceItems(p, m, pEnd, dictionary, transferSyntax, allowPrivate, recursionLevel + 1, tagList);
        end
    else
        if strcmp(vr, 'SQ')
            [pSequenceEnd, v] = parseSequenceItems(p, m, p + l - 1, dictionary, transferSyntax, allowPrivate, recursionLevel + 1, tagList);
            assert(pSequenceEnd == p + l);
        else
            v = m(p + (0:l-1));
        end        
        p = p + l;
    end
    
    %Ugly hack
    if strcmp(vr, 'UI')
        iNull = find(v == 0, 1);
        if ~isempty(iNull)
            v = v(1:iNull-1, 1);
        end
    end
    %End of ugly hack
    
    dataElement.t = tag;
    dataElement.vr = vr;
    dataElement.l = l;
    dataElement.v = v;
    dataElement.pStart = pStart;
    dataElement.pNextStart = p;
else
    if (tag == dicomHexPairToNumber('(FFFE,E00D)'))
        p = p + 8;
    end
    dataElement = [];
end

function [p, x] = parseSequenceItems(p, m, pEnd, dictionary, transferSyntax, allowPrivate, recursionLevel, tagList)

itemList = struct('t',{},'vr',{},'l',{},'v',{});

x = 1;

while ~isempty(x) && p <= pEnd
    [p, x] = parseSequenceItem(p, m, pEnd, dictionary, transferSyntax, allowPrivate, recursionLevel, tagList);
    if ~isempty(x)
        itemList(end + 1) = x;
    end
end

x = itemList;


function [p, item] = parseSequenceItem(p, m, pEnd, dictionary, transferSyntax, allowPrivate, recursionLevel, tagList)

padding = repmat(' ',[1 recursionLevel]);

groupNumber = typecast(m(p + (0:1)),'uint16');
elementNumber = typecast(m(p + (2:3)),'uint16');

%tag = sprintf('(%0.4X,%0.4X)', groupNumber, elementNumber);
tag = double(typecast([groupNumber, elementNumber],'uint32'));

% fprintf('%s%s %8s', padding, tag, 'ITEM');

assert(any(tag == cellfun(@dicomHexPairToNumber,{'(FFFE,E000)','(FFFE,E0DD)'})));

l = double(typecast(m(p + (4:7)),'uint32'));

assert(p + l - 1 <= pEnd || strcmp(sprintf('%.0x',l),'ffffffff'));

switch tag
    case dicomHexPairToNumber('(FFFE,E000)')
%         fprintf(' %8.0f:\n', l);
        
        %x = m(p + (0:l-1));
        [p, elementList] = parseDataElementsEarlyTermination(p + 8, m, p + l - 1, dictionary, transferSyntax, allowPrivate, recursionLevel + 1, tagList);
        
        item.t = tag;
        item.vr = 'ITEM';
        item.l = l;
        item.v = elementList;
    case dicomHexPairToNumber('(FFFE,E0DD)')
        p = p + 8;
        item = [];
end




