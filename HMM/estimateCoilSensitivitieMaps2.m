function [noiseWhiteningTransform, channelSensitivityMaps, adcChannelIds, debugStructure] = estimateCoilSensitivitieMaps2(s)

nAcquisitions = numel(unique(s.m{1}.mdh.sLC.ushAcquisition));

assert(nAcquisitions == 2, 'Expected two ICE acquisitions in calibration scan.');

assert(s.m{1}.mdh.ushSamplesInScan(end) == 0, 'Calibration scan terminated unexpectedly (incomplete data).');

nSets = numel(unique(s.m{1}.mdh.sLC.ushSet(s.m{1}.mdh.sLC.ushAcquisition == 0)));

assert(nSets == 2, 'Expected two ICE sets in calibration scan.');

ii1 = find(s.m{1}.mdh.sLC.ushAcquisition == 0 & s.m{1}.mdh.sLC.ushSet == 0 & s.m{1}.mdh.ushSamplesInScan == 128); % Surface coil calibration scan
ii2 = find(s.m{1}.mdh.sLC.ushAcquisition == 0 & s.m{1}.mdh.sLC.ushSet == 1 & s.m{1}.mdh.ushSamplesInScan == 128); % Body coil calibration scan
ii3 = find(s.m{1}.mdh.sLC.ushAcquisition == 0 & s.m{1}.mdh.sLC.ushSet == 0 & s.m{1}.mdh.ushSamplesInScan == 512); % First noise measurement.
ii4 = find(s.m{1}.mdh.sLC.ushAcquisition == 1 & s.m{1}.mdh.sLC.ushSet == 0 & s.m{1}.mdh.ushSamplesInScan == 512); % Second noise measurement.
ii5 = find(s.m{1}.mdh.ushSamplesInScan == 0); % End of scan packet.

n1 = numel(ii1);
n2 = numel(ii2);
n3 = numel(ii3);
n4 = numel(ii4);
n5 = numel(ii5);

assert(n1 == n2);
assert(n3 == n4);

nReadouts = size(s.m{1}.mdh,1);

assert(n1 + n2 + n3 + n4 + n5 == nReadouts);

x1 = s.m{1}.mdh.sLC.ushLine(ii1);
y1 = s.m{1}.mdh.sLC.ushPartition(ii1);

x2 = s.m{1}.mdh.sLC.ushLine(ii2);
y2 = s.m{1}.mdh.sLC.ushPartition(ii2);

% x3 = s.m{1}.mdh.sLC.ushLine(ii3);
% y3 = s.m{1}.mdh.sLC.ushPartition(ii3);
% 
% x4 = s.m{1}.mdh.sLC.ushLine(ii4);
% y4 = s.m{1}.mdh.sLC.ushPartition(ii4);

nSurfaceCoils = size(s.m{1}.mdh.channelHeaders{ii1(1)},1);
nBodyCoils = size(s.m{1}.mdh.channelHeaders{ii2(1)},1);

jj1 = sub2ind([32 32], x1 + 1, y1 + 1);
c1 = cell(32, 32);
c1(jj1) = cellfun(@(x)permute(cell2mat(x.rawdata'),[3 1 4 2]),s.m{1}.mdh.channelHeaders(ii1),'UniformOutput',false);
m1 = ~cellfun(@(x)isempty(x),c1);
c1(~m1) = {zeros(1,128,1,nSurfaceCoils,'single')};
C1 = cell2mat(permute(c1, [1 3 2 4]));

jj2 = sub2ind([32 32], x2 + 1, y2 + 1);
c2 = cell(32, 32);
c2(jj2) = cellfun(@(x)permute(cell2mat(x.rawdata'),[3 1 4 2]),s.m{1}.mdh.channelHeaders(ii2),'UniformOutput',false);
m2 = ~cellfun(@(x)isempty(x),c2);
c2(~m2) = {zeros(1,128,1,nBodyCoils,'single')};
C2 = cell2mat(permute(c2, [1 3 2 4]));

c3 = cell2mat(cellfun(@(x)permute(cell2mat(x.rawdata'),[3 1 4 2]),s.m{1}.mdh.channelHeaders(ii3),'UniformOutput',false));
C3 = reshape(c3, [128*512 nSurfaceCoils]);

c4 = cell2mat(cellfun(@(x)permute(cell2mat(x.rawdata'),[3 1 4 2]),s.m{1}.mdh.channelHeaders(ii4),'UniformOutput',false));
C4 = reshape(c4, [128*512 nSurfaceCoils]);

C34 = cat(1, C3, C4);

noiseCovarianceMatrix = C34'*C34/size(C34,1);

noiseColoringTransform = chol(noiseCovarianceMatrix);

noiseWhiteningTransform = inv(noiseColoringTransform);

C1decorr = tprod(C1, [1 2 3 -5], noiseWhiteningTransform, [-5 4]);

clear C1;

sigma1 = 15/2; % Cycles per FOV

[I, J, K] = ndgrid(-16:15,(-64:63)/2,-16:15);

F1 = exp(-(I.^2 + J.^2 + K.^2)/(2*sigma1^2));

baseResolution = double(s.m{1}.param.Meas{1}.MEAS.sKSpace.lBaseResolution);

nominalMapArraySize = [baseResolution, baseResolution, baseResolution];
reconstructionMapArraySize = [baseResolution, 2*baseResolution, baseResolution];

newArraySize1 = [2*reconstructionMapArraySize, size(C1decorr,4)];
newArraySize2 = [2*reconstructionMapArraySize, size(C2,4)];

sigma2 = 15/2;

[I, J, K] = ndgrid(-(newArraySize1(1)/2):(newArraySize1(1)/2-1),(-(newArraySize1(2)/2):(newArraySize1(2)/2-1))/2,-(newArraySize1(3)/2):(newArraySize1(3)/2-1));

F2 = exp(-(I.^2 + J.^2 + K.^2)/(2*sigma2^2));

D1ref = centeredFFT(centeredResize(bsxfun(@times, C1decorr, 1), newArraySize1), 3);
D2ref = centeredFFT(centeredResize(bsxfun(@times, C2, 1), newArraySize2), 3);

D1 = centeredFFT(centeredResize(bsxfun(@times, C1decorr, F1), newArraySize1), 3);
D1 = centeredFFT(centeredResize(bsxfun(@times, C1decorr, 1), newArraySize1), 3);

D2 = centeredFFT(centeredResize(bsxfun(@times, C2, F1), newArraySize2), 3);

% D1 = centeredFFT(C1, 3);
% D2 = centeredFFT(C2, 3);

%D1decorrref = tprod(D1ref, [1 2 3 -5], noiseWhiteningTransform, [-5 4]);

%clear D1;

D20 = sum(D2, 4); % I'm guessing here. Don't know what the two different BC channels mean. But it looks more homogenous if I sum them.
D20ref = sum(D2ref, 4);

J0 = bsxfun(@times, conj(D1), abs(D20));
W0 = sum(abs(D1).^2, 4);

JF0 = centeredFFT(bsxfun(@times, centeredIFFT(J0, 3), F2), 3);
WF0 = centeredFFT(bsxfun(@times, centeredIFFT(W0, 3), F2), 3);

% JF0 = gaussFilterArray(J0, [2 2 2]);
% WF0 = gaussFilterArray(W0, [2 2 2]);

C = bsxfun(@rdivide, JF0, WF0);

T2 = sum(C.*D1ref,4);

adcChannelIds = double(s.m{1}.mdh.channelHeaders{ii1(1)}.ushChannelId);

channelSensitivityMaps = ArrayVolume(C);

nominalFieldOfViewPRS = [s.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dPhaseFOV, s.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dReadoutFOV, s.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dThickness];

voxelSize = (nominalFieldOfViewPRS./nominalMapArraySize/2)';

channelSensitivityMaps.R = [0 0 1; -1 0 0; 0 -1 0]; % Ad hoc solution. 2016-08-26
channelSensitivityMaps.v = voxelSize;
channelSensitivityMaps.r0 = channelSensitivityMaps.R*(diag(voxelSize)*(-[floor(2*reconstructionMapArraySize'/2) + [0 0 -1]']));

debugStructure.D20ref = channelSensitivityMaps.cloneWithNewContent(D20ref);
debugStructure.T2 = channelSensitivityMaps.cloneWithNewContent(T2);
debugStructure.bodyCoilImage = channelSensitivityMaps.cloneWithNewContent(abs(D20ref));

% sigma = 15/1; % Cycles per FOV
% 
% [I, J, K] = ndgrid(-16:15,(-64:63)/2,-16:15);
% 
% F = exp(-(I.^2 + J.^2 + K.^2)/(2*sigma^2));
% 
% baseResolution = double(s.m{1}.param.Meas{1}.MEAS.sKSpace.lBaseResolution);
% 
% nominalMapArraySize = [baseResolution, baseResolution, baseResolution];
% reconstructionMapArraySize = [baseResolution, 2*baseResolution, baseResolution];
% 
% newArraySize1 = [2*reconstructionMapArraySize, size(C1,4)];
% newArraySize2 = [2*reconstructionMapArraySize, size(C2,4)];
% 
% D1 = centeredFFT(centeredResize(bsxfun(@times, C1, F), newArraySize1), 3);
% D2 = centeredFFT(centeredResize(bsxfun(@times, C2, F), newArraySize2), 3);
% 
% D1decorr = tprod(D1, [1 2 3 -5], noiseWhiteningTransform, [-5 4]);
% 
% D10 = sqrt(sum(abs(D1decorr).^2, 4));
% 
% D20 = sqrt(sum(abs(D2).^2, 4));
% 
% S1 = bsxfun(@rdivide, D20, D10);
% 
% finalMapArraySize = 2*nominalMapArraySize;
% 
% S12 = centeredResize(S1, finalMapArraySize);
% 
% D202 = centeredResize(D20, finalMapArraySize);
% D102 = centeredResize(D10, finalMapArraySize);
% 
% %S120 = sqrt(sum(abs(S12).^2,4));
% 
% a = ArrayVolume(S12);
% 
% fieldOfViewPRS = [s.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dPhaseFOV, s.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dReadoutFOV, s.m{1}.param.Meas{1}.MEAS.sSliceArray.asSlice{1}.dThickness];
% 
% voxelSize = (fieldOfViewPRS./finalMapArraySize)';
% 
% % R_physical_to_patient = convertVB17PatientPositionToGradientOrientation(s.m{1}.param.Meas{1}.YAPS.tPatientPosition);
% % 
% % R_logical_to_physical = quaternion_to_rotationmatrix(s.m{1}.mdh.sSD(1,:).aflQuaternion(1,:));
% % 
% % R_logical_to_patient = R_physical_to_patient*R_logical_to_physical;
% 
% a.R = [0 0 1; -1 0 0; 0 -1 0]; % Ad hoc solution. 2016-08-26
% a.v = voxelSize;
% a.r0 = a.R*(diag(voxelSize)*(-[floor(finalMapArraySize'/2) + [0 0 -1]']));
% 
% d202 = a.cloneWithNewContent(D202);
% d102 = a.cloneWithNewContent(D102);
% 
% coilCalibrationImages = {a, d202, d102};


























