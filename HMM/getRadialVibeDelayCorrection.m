function delayCorrection = getRadialVibeDelayCorrection(Z0, iiDelayCorrectionLines)

nReadoutSamples = size(Z0,1);

nDelayCorrectionLines = numel(iiDelayCorrectionLines);

assert(size(Z0,3) == nDelayCorrectionLines);

xi = ((1:size(Z0,1)) - (floor(size(Z0,1)/2) + 1))';

phi = mod(double(iiDelayCorrectionLines - 1)*(sqrt(5) - 1)/2*pi, 2*pi);

[~, ii] = sort(phi);

%phi2 = phi(ii);

Z2 = Z0(:,:,ii);

Q21 = fftshift(ifft(ifftshift(Z2,1),[],1),1);
Q22 = circshift((fftshift(ifft(ifftshift(conj(Z2),1),[],1),1)), nDelayCorrectionLines/2, 3);

Q2 = Q22.*Q21;

Q3 = squeeze(sum(Q2, 2));

Q3s = double(sum(Q3,2));

%thetaDiffSampled = angle(Q3s);
%c0 = 0.0011;
c0 = fmincon(@(x)sum(abs(Q3s - abs(Q3s).*exp(2*pi*1i*x*xi)).^2)/sum(abs(Q3s).^2),1/nReadoutSamples,[],[],[],[],-5, 5);

thetaDiffFitted = 2*pi*c0*xi;

%figure(178);
%plot(xi, angle(Q3s), xi, thetaDiffFitted);
%drawnow;

delayCorrection = exp(-1i*thetaDiffFitted/2);

%Q4 = bsxfun(@times, exp(-1i*thetaDiffFitted), Q3);

%Q4s = double(sum(Q4, 2));

%Z3 = fftshift(ifft(ifftshift(bsxfun(@times, exp(-1i*thetaDiffFitted/2), Q21),1),[],1),1);