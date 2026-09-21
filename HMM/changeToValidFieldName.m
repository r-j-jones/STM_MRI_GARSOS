function validFieldName = changeToValidFieldName(name, tag)

if isvarname(name)
    validFieldName = name;
else
    name(~(isstrprop(name,'alphanum') | name == '_')) = '_';
    name = name(1:min(63:end));
    if numel(name) == 0 || isstrprop(name(1),'digit')
        validFieldName = strcat(tag,name);
    end
end