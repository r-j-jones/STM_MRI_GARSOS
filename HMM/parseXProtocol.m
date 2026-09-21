function [x, t] = parseXProtocol(s)

x = readNodes(s, 1);

function [x, p] = readNodes(s, p)
    
p = skipWhitespace(s, p);

x = cell(1, 0);

moreNodes = true;

while moreNodes && p < numel(s)
    
    [y, pTmp] = readNode(s, p); p = pTmp;
    
    p = skipWhitespace(s, p);
    
    if ~isempty(y) || iscell(y) || ischar(y)
        x{end + 1} = y;
    else
        moreNodes = false;
    end
    
end

function [x, p2] = readNode(s, p1)

switch s(p1)
    case '<'
        [x, p2] = readTag(s, p1);
    case '{'
        [x, p2] = readBlock(s, p1);
    otherwise
        [x, p2] = readString(s, p1);
end

function [x, p] = readTag(s, p)

x = struct('tag','','name','','value',[],'characterNumber',0);

if s(p) == '<'
    x.characterNumber = p;
    p = p + 1;
    [x.tag, p] = readUntil(s, p, '.>');
    if s(p) == '.'
        p = p + 1;
        [x.name, p] = readUntil(s, p, '>');
        if(x.name(1) == '"' && x.name(1) == '"')
            x.name = x.name(2:end-1);
        end
    end
    p = p + 1;
    
    p = skipWhitespace(s, p);
    
    [x.value, p] = readNode(s, p);
    
    if ~iscell(x.value)
        x.value = {x.value};
    end
    
else
    x = [];
end

function [x, p] = readBlock(s, p)

% if p == 7938
%     1+1;
% end

if s(p) == '{'
    [x, p] = readNodes(s, p + 1);
else
    x = [];
end

if s(p) == '}'
    p = p + 1;
else
    error('Syntax error');
end

function [x, p] = readString(s, p)

  
switch s(p)
    case '"'        
        [x, p] = readUntil(s, p + 1, '"');        
        p = p + 1;
        while s(p) == '"'
            [xp, p] = readUntil(s, p + 1, '"');
            x = [x '"' xp]; %#ok
            p = p + 1;
        end
    case '}'
        x = [];
    otherwise
        if numel(s) >= p+16 && strcmp(s(p:p+16),'### ASCCONV BEGIN')
            [x, p] = readAsciiSerialization(s, p);
        else
            [x, pTmp] = readUntil(s, p, sprintf(' \t')); p = pTmp;
        end
end

function [x, p] = readAsciiSerialization(s, p)

assert(strcmp(s(p:p+16),'### ASCCONV BEGIN'),'ASCII serialization must begin with ### ASCCONV BEGIN')

k = strfind(s, '### ASCCONV END ###');

k = min(k(k > p));

p0 = p;

p = k + 19;

x = s(p0:p-1);

function p = skipWhitespace(s, p)

%ws = sprintf(' \t\r\n');
ws = char([32, 9, 13, 10]);
ns = numel(s);

if p <= ns
    while any(s(p) == ws)
        p = p + 1;
        if p > ns
            break;
        end
    end
end

function [ss, p] = readUntil(s, p, t)

p0 = p;

while p <= numel(s) && ~any(s(p) == t)
    p = p + 1;
end

ss = s(p0:p-1);












