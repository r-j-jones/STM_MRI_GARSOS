function channelHeaderFields = getChannelHeaderStructureVD11()

channelHeaderFields(1).name = 'ulTypeAndChannelLength';
channelHeaderFields(1).type = 'uint32';
channelHeaderFields(1).length = 4;

channelHeaderFields(2).name = 'lMeasUID';
channelHeaderFields(2).type = 'int32';
channelHeaderFields(2).length = 4;

channelHeaderFields(3).name = 'ulScanCounter';
channelHeaderFields(3).type = 'uint32';
channelHeaderFields(3).length = 4;

channelHeaderFields(4).name = 'ulReserved1';
channelHeaderFields(4).type = 'uint32';
channelHeaderFields(4).length = 4;

channelHeaderFields(5).name = 'ulSequenceTime';
channelHeaderFields(5).type = 'uint32';
channelHeaderFields(5).length = 4;

channelHeaderFields(6).name = 'ulUnused2';
channelHeaderFields(6).type = 'uint32';
channelHeaderFields(6).length = 4;

channelHeaderFields(7).name = 'ushChannelId';
channelHeaderFields(7).type = 'uint16';
channelHeaderFields(7).length = 2;

channelHeaderFields(8).name = 'ushUnused3';
channelHeaderFields(8).type = 'uint16';
channelHeaderFields(8).length = 2;

channelHeaderFields(9).name = 'ulCRC';
channelHeaderFields(9).type = 'uint32';
channelHeaderFields(9).length = 4;