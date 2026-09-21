function [manifest, manifestFile] = radial_vibe_18_5_4_timestamp_imFIAT(kSpaceID, scanDir,reconstructionDir, forceUpdate)

dstRootDir = reconstructionDir;

prefixString = 'MultiSeq';

dstSeriesDirName = sprintf('%s_Meta', prefixString);

dstDir = fullfile(dstRootDir, dstSeriesDirName);

permissiveMakeDirIfNotExist(dstDir);

manifestFile = fullfile(dstDir, 'deid-manifest.mat');

reconstructionTimestamp = '2017-02-20 16:06:50 -05:00';

[manifestLoaded, manifest, ~] = loadIfExistAndTimestampNewerThan(manifestFile, reconstructionTimestamp);

if ~manifestLoaded || forceUpdate
    
    manifest = struct();
    
    manifest.prefixString = prefixString;
    manifest.dstRootDir = dstRootDir;
    
    %rawDataFiles = examination.getRawDataFiles();
    rawDataDir = fullfile(scanDir, 'raw');
    rawDataFiles = dir(fullfile(rawDataDir, '*.dat'));      
    rawDataFiles = arrayfun(@(x)fullfile(rawDataDir, x.name), rawDataFiles, 'UniformOutput', false);
    
    % create a 6*1 cell, each file has one
    ss = cellfun(@(x)readVD11MultiRaidFileStructure(x, 'HeaderOnly', true), rawDataFiles, 'UniformOutput', false);

    sequenceIds = cellfun(@(x)x.m{end}.param.Meas{1}.HEADER.MeasUID, ss);

    selectedSequenceIds = [kSpaceID];
    assert(all(ismember(selectedSequenceIds, sequenceIds)));
    
    nSequences = numel(selectedSequenceIds);
    
    for iSequence = 1:nSequences
        
        iDceProtocol = find(sequenceIds == selectedSequenceIds(iSequence));
        
        %s = readVD11MultiRaidFileStructure(rawDataFiles{iDceProtocol},'CalibrationOnly',true);
        
        q = mapVBVD(rawDataFiles{iDceProtocol});
        
        nSpokes = double(q{end}.image.NLin * q{end}.image.NRep);
        
        % 0-3499 index
        iiSpokes = double(q{end}.image.Lin(q{end}.image.Par == q{end}.image.centerPar) - 1);
        % timestamp
        tSpokes = 2.5*double(q{end}.image.timestamp(q{end}.image.Par == q{end}.image.centerPar));
        
        % resort it
        [iiSpokes, iiOrder] = sort(iiSpokes);
        tSpokes = tSpokes(iiOrder);
        
        manifest.sequences(iSequence).sequenceId = sequenceIds(iDceProtocol);
        manifest.sequences(iSequence).rawDataFile = rawDataFiles{iDceProtocol};
        manifest.sequences(iSequence).nSpokes = nSpokes;
        manifest.sequences(iSequence).iiSpokes = iiSpokes;
        manifest.sequences(iSequence).spokeTimestamps = tSpokes;
        
    end
    
    saveWithTimestamp(manifestFile , manifest, iso8601Now());
    
    % if is already there, read the data from already has
    [manifestLoaded, manifest, ~] = loadIfExistAndTimestampNewerThan(manifestFile, reconstructionTimestamp);
    
    assert(manifestLoaded);
    
end


