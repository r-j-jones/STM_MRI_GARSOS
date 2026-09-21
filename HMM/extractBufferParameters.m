function p = extractBufferParameters(x)

switch class(x)
    case 'struct'
        p = extractXProtocolParameter(x);
    case 'char'
        p = extractSerializedAsciiParamters(x);
end

