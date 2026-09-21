function mdhFields = getMdhStructureVD11()

mdhFields(1).name = 'ulFlagsAndDMALength';
mdhFields(1).type = 'uint32';
mdhFields(1).length = 4;

mdhFields(2).name = 'lMeasUID';
mdhFields(2).type = 'int32';
mdhFields(2).length = 4;

mdhFields(3).name = 'ulScanCounter';
mdhFields(3).type = 'uint32';
mdhFields(3).length = 4;

mdhFields(4).name = 'ulTimeStamp';
mdhFields(4).type = 'uint32';
mdhFields(4).length = 4;

mdhFields(5).name = 'ulPMUTimeStamp';
mdhFields(5).type = 'uint32';
mdhFields(5).length = 4;

%Begin fix later
mdhFields(6).name = 'ushSystemType';
mdhFields(6).type = 'uint16';
mdhFields(6).length = 2;

mdhFields(7).name = 'ushPTABPosDelay';
mdhFields(7).type = 'uint16';
mdhFields(7).length = 2;

mdhFields(8).name = 'lPATBPosX';
mdhFields(8).type = 'int32';
mdhFields(8).length = 4;

mdhFields(9).name = 'lPATBPosY';
mdhFields(9).type = 'int32';
mdhFields(9).length = 4;

mdhFields(10).name = 'lPATBPosZ';
mdhFields(10).type = 'int32';
mdhFields(10).length = 4;

mdhFields(11).name = 'ulReserved1';
mdhFields(11).type = 'uint32';
mdhFields(11).length = 4;
%End fix later

mdhFields(6 + 6).name = 'aulEvalInfoMask';
mdhFields(6 + 6).type = 'uint32';
mdhFields(6 + 6).length = 8;

mdhFields(7 + 6).name = 'ushSamplesInScan';
mdhFields(7 + 6).type = 'uint16';
mdhFields(7 + 6).length = 2;

mdhFields(8 + 6).name = 'ushUsedChannels';
mdhFields(8 + 6).type = 'uint16';
mdhFields(8 + 6).length = 2;

mdhFields(9 + 6).name = 'sLC';
mdhFields(9 + 6).type(1).name = 'ushLine';
mdhFields(9 + 6).type(1).type = 'uint16';
mdhFields(9 + 6).type(1).length = 2;
mdhFields(9 + 6).type(2).name = 'ushAcquisition';
mdhFields(9 + 6).type(2).type = 'uint16';
mdhFields(9 + 6).type(2).length = 2;
mdhFields(9 + 6).type(3).name = 'ushSlice';
mdhFields(9 + 6).type(3).type = 'uint16';
mdhFields(9 + 6).type(3).length = 2;
mdhFields(9 + 6).type(4).name = 'ushPartition';
mdhFields(9 + 6).type(4).type = 'uint16';
mdhFields(9 + 6).type(4).length = 2;
mdhFields(9 + 6).type(5).name = 'ushEcho';
mdhFields(9 + 6).type(5).type = 'uint16';
mdhFields(9 + 6).type(5).length = 2;
mdhFields(9 + 6).type(6).name = 'ushPhase';
mdhFields(9 + 6).type(6).type = 'uint16';
mdhFields(9 + 6).type(6).length = 2;
mdhFields(9 + 6).type(7).name = 'ushRepetition';
mdhFields(9 + 6).type(7).type = 'uint16';
mdhFields(9 + 6).type(7).length = 2;
mdhFields(9 + 6).type(8).name = 'ushSet';
mdhFields(9 + 6).type(8).type = 'uint16';
mdhFields(9 + 6).type(8).length = 2;
mdhFields(9 + 6).type(9).name = 'ushSeg';
mdhFields(9 + 6).type(9).type = 'uint16';
mdhFields(9 + 6).type(9).length = 2;
mdhFields(9 + 6).type(10).name = 'ushIda';
mdhFields(9 + 6).type(10).type = 'uint16';
mdhFields(9 + 6).type(10).length = 2;
mdhFields(9 + 6).type(11).name = 'ushIdb';
mdhFields(9 + 6).type(11).type = 'uint16';
mdhFields(9 + 6).type(11).length = 2;
mdhFields(9 + 6).type(12).name = 'ushIdc';
mdhFields(9 + 6).type(12).type = 'uint16';
mdhFields(9 + 6).type(12).length = 2;
mdhFields(9 + 6).type(13).name = 'ushIdd';
mdhFields(9 + 6).type(13).type = 'uint16';
mdhFields(9 + 6).type(13).length = 2;
mdhFields(9 + 6).type(14).name = 'ushIde';
mdhFields(9 + 6).type(14).type = 'uint16';
mdhFields(9 + 6).type(14).length = 2;
mdhFields(9 + 6).length = 28;

mdhFields(10 + 6).name = 'sCutOff';
mdhFields(10 + 6).type(1).name = 'ushPre';
mdhFields(10 + 6).type(1).type = 'uint16';
mdhFields(10 + 6).type(1).length = 2;
mdhFields(10 + 6).type(2).name = 'ushPost';
mdhFields(10 + 6).type(2).type = 'uint16';
mdhFields(10 + 6).type(2).length = 2;
mdhFields(10 + 6).length = 4;

mdhFields(11 + 6).name = 'ushKSpaceCentreColumn';
mdhFields(11 + 6).type = 'uint16';
mdhFields(11 + 6).length = 2;

mdhFields(12 + 6).name = 'ushCoilSelect';
mdhFields(12 + 6).type = 'uint16';
mdhFields(12 + 6).length = 2;

mdhFields(13 + 6).name = 'fReadOutOffCentre';
mdhFields(13 + 6).type = 'single';
mdhFields(13 + 6).length = 4;

mdhFields(14 + 6).name = 'ulTimeSinceLastRF';
mdhFields(14 + 6).type = 'uint32';
mdhFields(14 + 6).length = 4;

mdhFields(15 + 6).name = 'ushKSpaceCentreLineNo';
mdhFields(15 + 6).type = 'uint16';
mdhFields(15 + 6).length = 2;

mdhFields(16 + 6).name = 'ushKSpaceCentrePartitionNo';
mdhFields(16 + 6).type = 'uint16';
mdhFields(16 + 6).length = 2;

% Maybe fix? SODA?
mdhFields(23).name = 'sSD';
mdhFields(23).type(1).name = 'sSlicePosVec';
mdhFields(23).type(1).type(1).name = 'flSag';
mdhFields(23).type(1).type(1).type = 'single';
mdhFields(23).type(1).type(1).length = 4;
mdhFields(23).type(1).type(2).name = 'flCor';
mdhFields(23).type(1).type(2).type = 'single';
mdhFields(23).type(1).type(2).length = 4;
mdhFields(23).type(1).type(3).name = 'flTra';
mdhFields(23).type(1).type(3).type = 'single';
mdhFields(23).type(1).type(3).length = 4;
mdhFields(23).type(1).length = 12;
mdhFields(23).type(2).name = 'aflQuaternion';
mdhFields(23).type(2).type = 'single';
mdhFields(23).type(2).length = 16;
mdhFields(23).length = 28;

mdhFields(24).name = 'aushIceProgramPara';
mdhFields(24).type = 'uint16';
mdhFields(24).length = 48;

mdhFields(25).name = 'aushReservedPara';
mdhFields(25).type = 'uint16';
mdhFields(25).length = 8;

mdhFields(26).name = 'ushApplicationCounter';
mdhFields(26).type = 'uint16';
mdhFields(26).length = 2;

mdhFields(27).name = 'ushApplicationMask';
mdhFields(27).type = 'uint16';
mdhFields(27).length = 2;

mdhFields(28).name = 'ulCRC';
mdhFields(28).type = 'uint32';
mdhFields(28).length = 4;

