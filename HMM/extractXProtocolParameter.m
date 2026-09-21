function p = extractXProtocolParameter(x)

% parameterTags = {'XProtocol','ParamMap','ParamArray','ParamBool','ParamLong','ParamString','ParamDouble'};

assert(strcmp(x.tag,'XProtocol'),'Input must be an XProtocol structure');

p = extractParameterNodes(x);

f = fields(p);

assert(numel(f) <= 1);

p = p.(f{1});

function p = extractParameterNodes(n)

switch n.tag
    case {'XProtocol', 'ParamMap', 'PipeService'}
        p = getParameterMapValue(n);
    case 'ParamString'
        p = getParameterStringValue(n);
    case {'ParamDouble', 'ParamLong'}
        p = getParameterDoubleValue(n);
    case 'ParamBool'
        p = getParameterBoolValue(n);
    case 'ParamArray'
        p = getParameterArrayValues(n);
    otherwise
        p = [];
end

function p = getParameterMapValue(n)

nodeNames = cellfun(@(x)x.name, n.value, 'UniformOutput', false);
nodeTags = cellfun(@(x)x.tag, n.value, 'UniformOutput', false);
m = strcmp(nodeNames, '');
nodeNames(m) = arrayfun(@(x,y)sprintf('%s%.0f',x{1},y),nodeTags(m),find(m),'UniformOutput',false);
nodeNames = cellfun(@(x,y)changeToValidFieldName(x,y),nodeNames, nodeTags, 'UniformOutput', false);
nodeContent = cellfun(@(x)extractParameterNodes(x), n.value, 'UniformOutput', false);
c = cellfun(@(x)~isempty(x)|ischar(x),nodeContent);
nCells = sum(c);
s = reshape([nodeNames(c); nodeContent(c)],[2*nCells 1]);
p = struct();
for i = 1:2:numel(s)
    p.(s{i}(1:min(63, end))) = s{i+1};
end
% p = struct(s{:});

function v = getParameterStringValue(n)
v = getParamValues(n);

if isempty(v)
    v = getDefaultParameterValues(n);
end

switch numel(v)
    case 1
        v = v{1};
    otherwise
        v = '';
end

function v = getParameterBoolValue(n)
v = getParamValues(n);

if isempty(v)
    v = getDefaultParameterValues(n);
end

switch numel(v)
    case 1
        v = logical(eval(v{1}));
    otherwise
        v = false;
end

function v = getParameterDoubleValue(n)
v = getParamValues(n);

if isempty(v)
    v = getDefaultParameterValues(n);
end

switch numel(v)
    case 0
        v = 0;
    otherwise
        v = cellfun(@(x)sscanf(x,'%f'),v);
end

function v = getParamValues(n)
m = cellfun(@(x)ischar(x),n.value);
v = n.value(m);
if numel(v) == 0
    v = getDefaultParameterValues(n);
end

function v = getDefaultParameterValues(n)
m = cellfun(@(x)isstruct(x)&&strcmp(x.tag,'Default'),n.value);
d = [n.value{m}];
switch numel(d)
    case 0
        v = {};
    case 1
        v = d.value;
    otherwise
        error('Multiple default values');
end

function v = getParameterArrayValues(n)

defaultNode = getDefaultParameterValues(n);

defaultNode = defaultNode{1};

% defaultValues = extractParameterNodes(defaultNode);

arrayValueMask = cellfun(@(x)iscell(x),n.value);

arrayValues = n.value(arrayValueMask);

for iElement = 1:numel(arrayValues);
%     if isempty(arrayValues{iElement})
%         arrayValues{iElement} = defaultValues;
%     else
        templateNode = defaultNode;
%         switch templateNode.tag
%             case 'ParamMap'
%                 for iMember = 1:numel(templateNode.value)
                    filledNode = distributeIntoChildren(templateNode, arrayValues{iElement});
%                 end
%             otherwise
%                 filledNode = templateNode;
%                 filledNode.value = arrayValues{iElement};
%         end
        arrayValues{iElement} = extractParameterNodes(filledNode);
%     end
end

v = arrayValues;

function filledNode = distributeIntoChildren(templateNode, values)

if ischar(values)
    filledNode = {values};
end

filledNode = templateNode;

if isstruct(filledNode)
    k = getParameterNodeIndices(filledNode.value);
    for iMember = 1:numel(values)
        if iMember <= numel(k)
            assignIndex = k(iMember);
            filledNode.value{assignIndex} = distributeIntoChildren(filledNode.value{assignIndex}, values{iMember});
        else
            assignIndex = numel(filledNode.value) + 1;
            filledNode.value{assignIndex} = values{iMember};
        end
    end
else
    filledNode = values;
end

function [k, ms] = getParameterNodeIndices(values)

assert(iscell(values));

parameterTagNames = {'ParamMap','ParamString','ParamDouble', 'ParamLong','ParamBool','ParamArray'};

ms = cellfun(@(x)(~isstruct(x))||any(strcmp(x.tag,parameterTagNames)), values);

k = find(ms);









