function [subjNbStr, experimentNbStr, scanStr] = parseDeidInputPath( inputPath )

[tmp,experimentNbStr] = fileparts(inputPath);
[~,subjNbStr] = fileparts(tmp);
scanStr = sprintf('%s_%s',subjNbStr,experimentNbStr);
