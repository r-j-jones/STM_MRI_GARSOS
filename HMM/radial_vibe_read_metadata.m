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
        
        manifest.sequences(iSequence).sequenceId = sequenceIds(iDceProtocol);
        manifest.sequences(iSequence).rawDataFile = rawDataFiles{iDceProtocol};
        manifest.sequences(iSequence).nSpokes = nSpokes;
        manifest.sequences(iSequence).iiSpokes = iiSpokes;
        manifest.sequences(iSequence).spokeTimestamps = tSpokes;
        manifest.sequences(iSequence).sqzSize = double(q{end}.image.sqzSize);
        manifest.sequences(iSequence).sqzDims = q{end}.image.sqzDims;


        %% Get sequence / acquisition info + params

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


        
    end

    fprintf(' - Saving manifest file: \n %s \n',manifestFile);
    saveWithTimestamp(manifestFile , manifest, iso8601Now());
    [manifestLoaded, ~, ~] = loadIfExistAndTimestampNewerThan(manifestFile, reconstructionTimestamp);
    assert(manifestLoaded);
    
end

end



%% HELPER FUNCTIONS

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
    voxelSizeReconstructed = (nominalFieldOfViewPRS./reconstructionMapArraySize/2)';
    
    %[ Store parameters to output struct
    paramlist2 = {'nAcquisitions','nSets','nReadouts','nSurfaceCoils',...
        'nBodyCoils','baseResolution','nominalMapArraySize','reconstructionMapArraySize',...
        'adcChannelIds','nominalFieldOfViewPRS','voxelSize','voxelSizeReconstructed'};
    
    for pp=1:length(paramlist2)
        currp  = paramlist2{pp};
        % eval(['kinfo.(currp) = ' currp ';']);
        kinfo.(currp) = eval(currp);
    end


end



