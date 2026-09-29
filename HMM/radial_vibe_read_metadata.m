function [manifest, manifestFile] = radial_vibe_read_metadata(kSpaceID, scanDir, outDir, forceUpdate)
%RADIAL_VIBE_READ_METADATA load raw data headers and extract metadata 
%
% Modified from radial_vibe_18_5_4_timestamp_imFIAT.m
% RJJ  |  09-25-2026

fprintf(' -- Running: radial_vibe_read_metadata -- \n');

if nargin<4
    forceUpdate = false;
end
if nargin<3
    error('Not enough input arguments');
end

dstRootDir = outDir;
prefixString = 'RadialVibe';
dstSeriesDirName = sprintf('%s_Meta', prefixString);
dstDir = fullfile(dstRootDir, dstSeriesDirName);
permissiveMakeDirIfNotExist(dstDir);

manifestFileName = 'garsos_meta-manifest.mat';
manifestFile = fullfile(dstDir, manifestFileName);
fprintf(' -Manifest file: %s \n',manifestFile);

reconstructionTimestamp = '2017-02-20 16:06:50 -05:00';
[manifestLoaded, manifest, ~] = loadIfExistAndTimestampNewerThan(manifestFile, reconstructionTimestamp);

fprintf(' -Manifest loaded = %g \n',manifestLoaded);
fprintf(' -Force update    = %g \n',forceUpdate);

if ~manifestLoaded || forceUpdate

    fprintf(' - Running main loop now \n');
    
    manifest = struct();
    
    manifest.prefixString = prefixString;
    manifest.dstRootDir = dstRootDir;
    manifest.manifestFileName = manifestFileName;
    
    rawDataDir = fullfile(scanDir, 'raw');
    rawDataFiles = dir(fullfile(rawDataDir, '*.dat'));      
    rawDataFiles = arrayfun(@(x)fullfile(rawDataDir, x.name), rawDataFiles, 'UniformOutput', false);
    
    ss = cellfun(@(x)readVD11MultiRaidFileStructure(x, 'HeaderOnly', true), rawDataFiles, 'UniformOutput', false);

    sequenceIds = cellfun(@(x)x.m{end}.param.Meas{1}.HEADER.MeasUID, ss);
    selectedSequenceIds = [kSpaceID]; %#ok<NBRAK2>
    assert(all(ismember(selectedSequenceIds, sequenceIds)));
    nSequences = numel(selectedSequenceIds);
    manifest.nSequences = nSequences;
    
    for iSequence = 1:nSequences

        fprintf(' - on sequence %d / %d \n',iSequence,nSequences);
        
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

    fprintf(' - Saving manifest file: \n %s \n',manifestFile);
    saveWithTimestamp(manifestFile , manifest, iso8601Now());
    [manifestLoaded, ~, ~] = loadIfExistAndTimestampNewerThan(manifestFile, reconstructionTimestamp);
    assert(manifestLoaded);
    
end


