function [kSpaceID,numseq] = parseKspaceIds( kSpaceID1, kSpaceID2 )
% function [kSpaceID,numseq] = parseKspaceIds( kSpaceID1, kSpaceID2 )

numseq = 2;
kSpaceID = [kSpaceID1, kSpaceID2];
if kSpaceID2 == -1
    kSpaceID = kSpaceID1;
    numseq = 1;        
end