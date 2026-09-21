function p = extractSerializedAsciiParamters(s) %#ok

r = strsplit(s, sprintf('\n'));

r = strtrim(r);

assert(strcmp(r{1}(1:17),'### ASCCONV BEGIN') & strcmp(r{end},'### ASCCONV END ###'), 'Not proper serialized ASCII.')

r = strcat('p.MEAS.',r(2:end-1),';');
r = strrep(r,'"','''');
r = regexprep(r,'\[(?<index>[0-9]+)\](?<tail>\S)','\{${sprintf(''%.0f'',str2num($<index>) + 1)}\}$<tail>');
r = regexprep(r,'\[(?<index>[0-9]+)\]\s','\(${sprintf(''%.0f'',str2num($<index>) + 1)}\)');
r = regexprep(r,'0x(?<number>[0-9a-f])+','${sprintf(''%.0f'',hex2dec($<number>))}');
r = regexprep(r,'.*__.*','');
%r = regexprep(r,'#.*;',';');

eval(strcat(r{:})); % The output parameter p comes from here.