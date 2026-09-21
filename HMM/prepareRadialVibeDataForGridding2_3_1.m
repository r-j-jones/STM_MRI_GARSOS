function [sampledData, kSpaceLocations, densityCompensation, iiSpoke, timeStamp, nPartitions, nLines, nSamples] = prepareRadialVibeDataForGridding2_3_1(q, s, noiseWhiteningTransform, channelIds)

%adcChannelIds = double(s.m{end}.mdh.channelHeaders{1}.ushChannelId);

%modified by Yuhang
% qq = squeeze(q.image());
%         [a,b,lin,c,rep]=size(qq);
%         q2 = qq(:,:,1,:,1);
%         for ii = 1:rep
%             for jj = 1:lin
%             q2(:,:,end+1,:) = squeeze(qq(:,:,jj,:,ii));
%             end
%         end
%         q2 = q2(:,:,2:end,:);
%         nAcq = q.image.NAcq;
%         nPar = q.image.NPar;
%         Lin2 = q.image.Lin;
%         for iAcq=1:nAcq
%             Lin2(iAcq)=ceil(iAcq/nPar);            
%         end
qq = squeeze(q.image());
        [a,b,lin,c,rep]=size(qq);
        q2 = zeros(a,b,lin*rep,c);
        for ii = 1:rep
            for jj = 1:lin
             q2(:,:,(ii-1)*lin+jj,:) = qq(:,:,jj,:,ii);
            end
        end
        q2=single(q2);
        nAcq = q.image.NAcq;
        nPar = q.image.NPar;
        Lin2 = q.image.Lin;
        for iAcq=1:nAcq
            Lin2(iAcq)=ceil(iAcq/nPar);            
        end
%%%%%%%

nTotalLines = double(q.image.NAcq);

nPreparationLines = 0;

nLines = double(rep*lin);
%nLines = q.image.NLin;
 
nCollectedPartitions = double(q.image.NPar);
nPartitions = s.m{end}.param.Meas{1}.MEAS.sKSpace.lPartitions;

assert(nTotalLines == nLines*nCollectedPartitions + nPreparationLines);

%validRows = m.mdh(nPreparationLines + 1:end-1,:);

nSamples = double(q.image.NCol);

%assert(numel(nSamples) == 1);

[iiPartition, iiLine] = ndgrid(0:nPartitions-1, 0:nLines-1);

%cPartition = floor(s.m{end}.param.Meas{1}.MEAS.sKSpace.lPartitions/2) - (s.m{end}.param.Meas{1}.MEAS.sKSpace.lPartitions - nPartitions);
cPartition = floor(s.m{end}.param.Meas{1}.MEAS.sKSpace.lPartitions/2);

rho = ((1:nSamples) - (floor(nSamples/2) + 1))/nSamples;
deltaRho = 1/nSamples;

%phi = double(validRows.sLC.ushLine)*pi*(3 - sqrt(5));
%phi = double(validRows.sLC.ushLine)*(sqrt(5) - 1)/2*pi;

%  phi = double(iiLine(:))*(sqrt(5) - 1)/2*pi;

  %iiLine2 = iiLine(:,18*lin+1:19*lin);
%  iiLine2 = iiLine(:,end-2*lin+1:end-lin);
  iiLine2 = iiLine(:,lin+1:2*lin);
  iiLine3 = repmat(iiLine2,1,rep);
  phi = double(iiLine3(:))*(sqrt(5) - 1)/2*pi;

  iiLine2b = iiLine(:,1:lin);
  iiLine3b = repmat(iiLine2b,1,rep);
  phib = double(iiLine3b(:))*(sqrt(5) - 1)/2*pi;

  

kx = cos(phi);
ky = sin(phi);
%kx = - cos(phi);
%ky = - sin(phi);

kx = kx*rho;
ky = ky*rho;
kz = (double(iiPartition(:)) - cPartition)/s.m{end}.param.Meas{1}.MEAS.sKSpace.lPartitions*s.m{end}.param.Meas{1}.MEAS.sKSpace.dSliceResolution;
%deltaZ = 1/s.m{end}.param.Meas{1}.MEAS.sKSpace.lPartitions*s.m{end}.param.Meas{1}.MEAS.sKSpace.dSliceResolution;
kz = repmat(kz, [1 nSamples]);

kx = reshape(transpose(kx),[nLines*nPartitions*nSamples 1]);
ky = reshape(transpose(ky),[nLines*nPartitions*nSamples 1]);
kz = reshape(transpose(kz),[nLines*nPartitions*nSamples 1]);

iiSpoke = repmat(double(iiLine(:)),[1 nSamples])';
iiSpoke = reshape(iiSpoke, [nLines*nPartitions*nSamples 1]);

% timeStamp = repmat(double(validRows.ulTimeStamp),[1 nSamples])';
% timeStamp = reshape(timeStamp, [nLines*nPartitions*nSamples 1]);

timeStamp = [];

%nChannels = q.image.NCha;

nFourierPartitions = double(nPartitions - q.image.NPar);

centerPartition = unique(q.image.centerPar);

nDelayCorrectionLines = floor(nLines/4)*2;

iiDelayCorrectionLines = nLines-nDelayCorrectionLines+1:nLines;

%Z0 = squeeze(q.image(:, :, nLines-nDelayCorrectionLines+1:nLines, centerPartition,:));
Z0 = squeeze(q2(:, :, nLines-nDelayCorrectionLines+1:nLines, centerPartition));

delayCorrection = getRadialVibeDelayCorrection(Z0, iiDelayCorrectionLines);

%delayCorrection = ones(size(delayCorrection));

%Z = q.image(:,:,:,:);
Z = q2;

Z = fftshift(fft(bsxfun(@times, ifft(ifftshift(Z,1),[],1), ifftshift(delayCorrection,1)),[],1),1);

Z = permute(Z, [1 4 3 2]);
Z = padarray(Z,[0 nFourierPartitions 0 0],0,'pre');
%Z = partialFourierWithPocs(Z, 2, 1:nFourierPartitions, [1 2]);
Z = decorrelateCoils2(Z, noiseWhiteningTransform); 
nTransformedChannels = size(noiseWhiteningTransform, 2);
Z = reshape(Z, [nLines*nPartitions*nSamples, nTransformedChannels]);

%sampledData = cell2mat(cellfun(@(x)cell2mat(x.rawdata'),validRows.channelHeaders,'UniformOutput', false));

sampledData = double(Z);

sampledData = double(sampledData);
kSpaceLocations = [kx ky kz];

RHOxy = sqrt(sum(kSpaceLocations(:,1:2).^2, 2));

upperRHOxy = (RHOxy + deltaRho/2);
lowerRHOxy = (RHOxy - deltaRho/2);

densityCompensation = pi*abs(upperRHOxy.^2 - sign(lowerRHOxy.*upperRHOxy).*lowerRHOxy.^2)/2; % Area of half a shell sector (other half for other side of spoke).

%densityCompensation = pi*abs(upperRHOxy.^2 - sign(lowerRHOxy.*upperRHOxy).*lowerRHOxy.^2); % Area of a whole shell sector (wrong!). Using anyway for consistency. All images twice as bright as they should be.

%densityCompensation = sqrt(sum(kSpaceLocations(:,1:2).^2, 2));





