function y = dicomHexPairToNumber(x)

y = double(typecast(uint16(sscanf(x,'(%x,%x)')),'uint32'));