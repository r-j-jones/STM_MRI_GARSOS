function [manifest, manifestFile] = radial_vibe_get_acq_params(kSpaceID, scanDir,reconstructionDir, forceUpdate)

dstRootDir = reconstructionDir;

prefixString = 'MultiSeq';

dstSeriesDirName = sprintf('%s_Meta', prefixString);

dstDir = fullfile(dstRootDir, dstSeriesDirName);

permissiveMakeDirIfNotExist(dstDir);

manifestFile = fullfile(dstDir, 'deid-AcqInfo-manifest.mat');

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
        
        s_calonly = readVD11MultiRaidFileStructure(rawDataFiles{iDceProtocol},'CalibrationOnly',true);
        
        q = mapVBVD(rawDataFiles{iDceProtocol});
        qhdr = q{end}.hdr;
        
        nSpokes = double(q{end}.image.NLin * q{end}.image.NRep);
        
        % 0-3499 index
        iiSpokes = double(q{end}.image.Lin(q{end}.image.Par == q{end}.image.centerPar) - 1);
        % timestamp
        tSpokes = 2.5*double(q{end}.image.timestamp(q{end}.image.Par == q{end}.image.centerPar));
        
        % resort it
        [iiSpokes, iiOrder] = sort(iiSpokes);
        tSpokes = tSpokes(iiOrder);
        
        % STORE
        manifest.sequences(iSequence).sequenceId = sequenceIds(iDceProtocol);
        manifest.sequences(iSequence).rawDataFile = rawDataFiles{iDceProtocol};
        manifest.sequences(iSequence).nSpokes = nSpokes;
        manifest.sequences(iSequence).iiSpokes = iiSpokes;
        manifest.sequences(iSequence).spokeTimestamps = tSpokes;
        manifest.sequences(iSequence).sqzSize = double(q{end}.image.sqzSize);
        manifest.sequences(iSequence).sqzDims = q{end}.image.sqzDims;


        ConfigFields = {'Interleaves','NAve','NAveMeas','NImageCols','NImageLins','NImagePar',...
            'NPar','NParMeas','NLin','NLinMeas',...
            'NTruePar','NoImagesPerSlab', 'NoOfFourierPartitions','NRoFtlen', 'PeFOV', ...
            'PhaseEncodingLines', 'PhaseFoV', 'PhasePartialFourierFactor', 'ReadFoV', ...
            'ReadoutOversamplingFactor', 'RoFOV', 'SequenceDescription', ...
            'TR', 'TimePerSpoke_us','PatientSex'};
        for ind=1:length(ConfigFields)
            currField = ConfigFields{ind};
            manifest.sequences(iSequence).hdr.(currField) = qhdr.Config.(currField);
        end

        MeasYapsFields = {'adFlipAngleDegree','alTR','alTE'};
        for ind=1:length(MeasYapsFields)
            currField = MeasYapsFields{ind};
            manifest.sequences(iSequence).hdr.(currField) = qhdr.MeasYaps.(currField){1};
        end

        DicomFields = {'flPatientAge','adFlipAngleDegree','flUsedPatientWeight'};
        for ind=1:length(DicomFields)
            currField = DicomFields{ind};
            manifest.sequences(iSequence).hdr.(currField) = qhdr.Dicom.(currField);
        end

        [ kinfo ] = returnKspaceImagingParameters(q{end}, ss{iDceProtocol}, s_calonly);
        manifest.sequences(iSequence).kinfo = kinfo;
        
        % manifest.sequences(iSequence).hdr.Interleaves = qhdr.Config.Interleaves;
        % manifest.sequences(iSequence).hdr.NAve = qhdr.Config.NAve;
        % manifest.sequences(iSequence).hdr.NAveMeas = qhdr.Config.NAveMeas;
        % manifest.sequences(iSequence).hdr.NImageCols = qhdr.Config.NImageCols;
        % manifest.sequences(iSequence).hdr.NImageLins = qhdr.Config.NImageLins;
        % manifest.sequences(iSequence).hdr.NImagePar = qhdr.Config.NImagePar;
        % manifest.sequences(iSequence).hdr.NPar = qhdr.Config.NPar;
        % manifest.sequences(iSequence).hdr.NParMeas = qhdr.Config.NParMeas;
        % manifest.sequences(iSequence).hdr.NLin = qhdr.Config.NLin;
        % manifest.sequences(iSequence).hdr.NLinMeas = qhdr.Config.NLinMeas;
        
    % Config 
        % 'NTruePar','NoImagesPerSlab', 'NoOfFourierPartitions','NRoFtlen', 'PeFOV', ...
        %     'PhaseEncodingLines', 'PhaseFoV', 'PhasePartialFourierFactor', 'ReadFoV', ...
        %     'ReadoutOversamplingFactor', 'RoFOV', 'SequenceDescription', ...
        %     'TR', 'TimePerSpoke_us'
    % MeasYaps
        % adFlipAngleDegree{1}
        % alTR{1}
        % alTE{1}
    % Dicom
        % flPatientAge  
        % adFlipAngleDegree
        % flUsedPatientWeight  
    end
    
    saveWithTimestamp(manifestFile , manifest, iso8601Now());
    
    % if is already there, read the data from already has
    [manifestLoaded, manifest, ~] = loadIfExistAndTimestampNewerThan(manifestFile, reconstructionTimestamp);
    
    assert(manifestLoaded);
    
end


end


function [ kinfo ] = returnKspaceImagingParameters(q_, s, s_calonly)


kinfo = [];


%% [ FROM prepareRadialVibeData...

%adcChannelIds = double(s.m{end}.mdh.channelHeaders{1}.ushChannelId);

nTotalLines = double(q_.image.NAcq);
nPreparationLines = 0;
nLines = double(q_.image.NLin);
nCollectedPartitions = double(q_.image.NPar);
nPartitions = s.m{end}.param.Meas{1}.MEAS.sKSpace.lPartitions;

phasePF = s.m{end}.param.Meas{1}.MEAS.sKSpace.ucPhasePartialFourier;
slicePF = s.m{end}.param.Meas{1}.MEAS.sKSpace.ucSlicePartialFourier;
readoutPF = s.m{end}.param.Meas{1}.MEAS.sKSpace.ucReadoutPartialFourier;

% assert(nTotalLines == nLines*nCollectedPartitions + nPreparationLines);
%validRows = m.mdh(nPreparationLines + 1:end-1,:);

nSamples = double(q_.image.NCol);
% rho = ((1:nSamples) - (floor(nSamples/2) + 1))/nSamples;
% deltaRho = 1/nSamples;

nChannels = q_.image.NCha;
nFourierPartitions = double(nPartitions - q_.image.NPar);
centerPartition = unique(q_.image.centerPar);
nDelayCorrectionLines = floor(q_.image.NLin/4)*2;

% nTransformedChannels = size(noiseWhiteningTransform, 2);

paramlist1 = {'nTotalLines','nPreparationLines','nLines','nCollectedPartitions',...
    'nPartitions','nSamples','nChannels','nFourierPartitions',...
    'centerPartition','nDelayCorrectionLines','phasePF','slicePF','readoutPF'}; %,'nTransformedChannels'};

for pp=1:length(paramlist1)
    currp  = paramlist1{pp};
    % eval(['kinfo.(currp) = ' currp ';']);
    kinfo.(currp) = eval(currp);
end


%% [ FROM ESTIMATE_COIL_SENSITIVITIES2

nAcquisitions = numel(unique(s_calonly.m{1}.mdh.sLC.ushAcquisition));

assert(nAcquisitions == 2, 'Expected two ICE acquisitions in calibration scan.');

assert(s_calonly.m{1}.mdh.ushSamplesInScan(end) == 0, 'Calibration scan terminated unexpectedly (incomplete data).');

nSets = numel(unique(s_calonly.m{1}.mdh.sLC.ushSet(s_calonly.m{1}.mdh.sLC.ushAcquisition == 0)));

assert(nSets == 2, 'Expected two ICE sets in calibration scan.');

ii1 = find(s_calonly.m{1}.mdh.sLC.ushAcquisition == 0 & s_calonly.m{1}.mdh.sLC.ushSet == 0 & s_calonly.m{1}.mdh.ushSamplesInScan == 128); % Surface coil calibration scan
ii2 = find(s_calonly.m{1}.mdh.sLC.ushAcquisition == 0 & s_calonly.m{1}.mdh.sLC.ushSet == 1 & s_calonly.m{1}.mdh.ushSamplesInScan == 128); % Body coil calibration scan
ii3 = find(s_calonly.m{1}.mdh.sLC.ushAcquisition == 0 & s_calonly.m{1}.mdh.sLC.ushSet == 0 & s_calonly.m{1}.mdh.ushSamplesInScan == 512); % First noise measurement.
ii4 = find(s_calonly.m{1}.mdh.sLC.ushAcquisition == 1 & s_calonly.m{1}.mdh.sLC.ushSet == 0 & s_calonly.m{1}.mdh.ushSamplesInScan == 512); % Second noise measurement.
ii5 = find(s_calonly.m{1}.mdh.ushSamplesInScan == 0); % End of scan packet.

n1 = numel(ii1);
n2 = numel(ii2);
n3 = numel(ii3);
n4 = numel(ii4);
n5 = numel(ii5);

assert(n1 == n2);
assert(n3 == n4);

nReadouts = size(s_calonly.m{1}.mdh,1);

assert(n1 + n2 + n3 + n4 + n5 == nReadouts);

x1 = s_calonly.m{1}.mdh.sLC.ushLine(ii1);
y1 = s_calonly.m{1}.mdh.sLC.ushPartition(ii1);

x2 = s_calonly.m{1}.mdh.sLC.ushLine(ii2);
y2 = s_calonly.m{1}.mdh.sLC.ushPartition(ii2);

% x3 = s_calonly.m{1}.mdh.sLC.ushLine(ii3);
% y3 = s_calonly.m{1}.mdh.sLC.ushPartition(ii3);
% 
% x4 = s_calonly.m{1}.mdh.sLC.ushLine(ii4);
% y4 = s_calonly.m{1}.mdh.sLC.ushPartition(ii4);

nSurfaceCoils = size(s_calonly.m{1}.mdh.channelHeaders{ii1(1)},1);
nBodyCoils = size(s_calonly.m{1}.mdh.channelHeaders{ii2(1)},1);
baseResolution = double(s_calonly.m{1}.param.Meas{1}.MEAS.sKSpace.lBaseResolution);
nominalMapArraySize = [baseResolution, baseResolution, baseResolution];
reconstructionMapArraySize = [baseResolution, 2*baseResolution, baseResolution];
adcChannelIds = double(s_calonly.m{1}.mdh.channelHeaders{ii1(1)}.ushChannelId);
nominalFieldOfViewPRS = [s_calonly.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dPhaseFOV, ...
    s_calonly.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dReadoutFOV, ...
    s_calonly.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dThickness];
voxelSize = (nominalFieldOfViewPRS./nominalMapArraySize/2)';
% voxelSizeReconstructed = (nominalFieldOfViewPRS./reconstructionMapArraySize/2)';

%[ Store parameters to output struct
paramlist2 = {'nAcquisitions','nSets','nReadouts','nSurfaceCoils',...
    'nBodyCoils','baseResolution','nominalMapArraySize','reconstructionMapArraySize',...
    'adcChannelIds','nominalFieldOfViewPRS'};

for pp=1:length(paramlist2)
    currp  = paramlist2{pp};
    % eval(['kinfo.(currp) = ' currp ';']);
    kinfo.(currp) = eval(currp);
end


end
