% PET reconstruct from list detector responses using Fourier rebin(FORE).
% Written by Shiwei Zhou, Youfang Lai 07/06/2020
% Reference: 1. Basics of PET Imaging Physics, Chemistry, and Regulations Second Edition
%               Page 43-46, 71-73
%            2. Defrise, Michel, et al. "Exact and approximate rebinning algorithms for 3-D PET data."
%               IEEE transactions on medical imaging 16.2 (1997): 145-158.
% Reference function: https://www.mathworks.com/help/images/ref/iradon.html

% Open bainary files, which have three columns, representing x, y and z
fid=fopen('coin5.dat','rb');
detectorPositions=fread(fid,'float64');
fclose(fid);
detectorPositions=reshape(detectorPositions,[],6);
% figure(1);
% clf;
% scatter(f1(:,1),f1(:,2));

deviceR = 80.5;
deviceNRing = 80;
deviceSigma = 1.59;
ddelta = deviceSigma/2/deviceR;
deviceL = 127;
MRD = 79;    
span = 3;

%     data(:,1) = data(:,1)-xOffset(kkk);
%     data(:,4) = data(:,4)-xOffset(kkk);
%     xOffset(kkk) = 0;    

% Change the dector positions to Line of Response.
rxy = data(:,1:3)-data(:,4:6);
tmp = rxy(:,3)~=0;
rxy(tmp,:) = rxy(tmp,:).*sign(rxy(tmp,3));
cosphi = rxy(:,1)./sqrt(sum(rxy(:,1:2).^2,2));
phi = acos(cosphi)*180/pi + 180*(rxy(:,2)<0);
k = rxy(:,2)./rxy(:,1);
r = (-k.*data(:,1)+data(:,2)).*cosphi;
r(isinf(k)) = data(isinf(k),1);
r(phi>180) = -r(phi>180);
phi(phi>180) = phi(phi>180)-180;
%     delta = rxy(:,3)./sqrt(sum(rxy(:,1:2).^2,2));
%     zmid = (data(:,3)+data(:,6))/2;

di = floor((data(:,3)+deviceL/2)/deviceSigma);
dj = floor((data(:,6)+deviceL/2)/deviceSigma);
zmid = -(deviceL-deviceSigma)/2+(di+dj)*deviceSigma/2;
dij = di -dj;
dunderscore = span*floor((dij+(span-1)/2)/span);
delta = dunderscore*ddelta;

detectorLines = [r,phi,zmid,delta];

% Construct the sinogram
rMin = -deviceR;
rMax = deviceR;
rSize = rMax - rMin;
%     rNumber = 2^10;
rNumber = 1660;
rRange = linspace(rMin , rMax ,rNumber);
omegaRange = 2 * pi * 1 / (rRange(2) - rRange(1)) * (-rNumber/2 : rNumber/2 -1) /rNumber;

phiNumber = 256;
phiRange = linspace(0,180, phiNumber);
kRange = 2 * pi * 1 / (phiRange(2) - phiRange(1)) * (-phiNumber/2 : phiNumber/2 -1) /phiNumber;

zDirectionMax = 24;
zDirectionMin = -24;    
zDirectionNumber = 61;
zDirectionRange = linspace(zDirectionMin , zDirectionMax ,zDirectionNumber);
zDirectionDifference = zDirectionRange(2) - zDirectionRange(1);

deltaMax = 0.36;%(MRD+1)/deviceNRing*deviceL/(rMax-rMin);
deltaMin = -deltaMax;
deltaNumber = ceil((deviceNRing+1)/span);
deltaRange = linspace(deltaMin , deltaMax ,deltaNumber);

sinogramImages = zeros(length(rRange), length(phiRange), length(zDirectionRange), length(deltaRange));
for idx = 1 : size(detectorLines,1)
    r = detectorLines(idx , 1);
    phi = detectorLines(idx , 2);
    z = detectorLines(idx , 3);
    delta = detectorLines(idx , 4);
    if ( r > rMax || r < rMin || z > zDirectionMax || z < zDirectionMin || delta > deltaMax || delta < deltaMin)
        continue;
    end
    rPosition = floor((r - rMin) / (rRange(2) - rRange(1))) +1;
    phiPosition = floor(phi/(phiRange(2)-phiRange(1)))+ 1;
    zPosition = floor((z - zDirectionMin) / zDirectionDifference) +1;
    deltaPosition = floor((delta - deltaMin) / (deltaRange(2)-deltaRange(1))) +1;
    sinogramImages(rPosition , phiPosition , zPosition , deltaPosition) = sinogramImages(rPosition , phiPosition , zPosition , deltaPosition) + 1;
end
clear idx r phi z delta rPosition phiPosition zPosition deltaPosition;

% Fourier rebin
deltaLim = span*ddelta;
kLim = 2;
rOmega = 70;
omegaLim = 2*(omegaRange(2) - omegaRange(1));
sinogramRebinFFTImages = zeros(length(rRange), length(phiRange), length(zDirectionRange));
weightForSino = zeros(length(rRange), length(phiRange), length(zDirectionRange));
for i1 = 1 : size(sinogramImages , 3)  % sweep all z
    if ( sum(abs(sinogramImages(: , : , i1 , :)),'all') ==0)       % for fast calc, skip this one
        continue;
    end
    for i2 = 1 : size(sinogramImages , 4) % sweep all delta
        if ( sum(abs(sinogramImages(: , : , i1 , i2)),'all') ==0)       % for fast calc, skip this one
            continue;
        end
        sinogramFFTImage = fftshift(fft2(sinogramImages(: , : , i1 , i2)));
        for i3 = 1 : size(sinogramImages , 1 )   % sweep all omega
            for i4 = 1 : size(sinogramImages , 2 )   % sweep all k
                delta1 = min([(deviceL/2 - zDirectionRange(i1)) / (deviceR + kRange(i4)/omegaRange(i3)) , ...
                               (deviceL/2 + zDirectionRange(i1)) / (deviceR - kRange(i4)/omegaRange(i3)) , ...
                               deviceL/2/deviceR]);                       
                delta2 = min([ (deviceL/2 - abs(zDirectionRange(i1)))/deviceR , deltaLim]);   
                if( abs(deltaRange(i2))<=abs(delta1) && abs(kRange(i4) / omegaRange(i3))<rOmega && (abs(kRange(i4))>kLim ||  abs(omegaRange(i3)) >omegaLim) ) % region 1
                    zShift = zDirectionRange(i1) - deltaRange(i2) * kRange(i4) / omegaRange(i3);
                    zShiftUpIndex = floor((zShift - zDirectionMin) / zDirectionDifference ) + 2;
                    zShiftDownIndex = floor((zShift - zDirectionMin) / zDirectionDifference ) + 1;                        
                    if(zShiftDownIndex<zDirectionNumber+1 &&zShiftDownIndex>0)
                        m = (zShift - zDirectionRange(zShiftDownIndex))/zDirectionDifference;
                        sinogramRebinFFTImages(i3 , i4 , zShiftDownIndex) = sinogramRebinFFTImages(i3 , i4 , zShiftDownIndex) + (1-m) * sinogramFFTImage(i3 , i4);
                        weightForSino(i3 , i4 , zShiftDownIndex) = weightForSino(i3 , i4 , zShiftDownIndex)+(1-m);
                    end
                    if(zShiftUpIndex<zDirectionNumber+1 && zShiftUpIndex>1)
                        m = (zShift - zDirectionRange(zShiftDownIndex))/zDirectionDifference;
                        sinogramRebinFFTImages(i3 , i4 , zShiftUpIndex) = sinogramRebinFFTImages(i3 , i4 , zShiftUpIndex) + m * sinogramFFTImage(i3 , i4);
                        weightForSino(i3 , i4 , zShiftDownIndex) = weightForSino(i3 , i4 , zShiftDownIndex)+m;
                    end
                end
                if (abs(deltaRange(i2))<=abs(delta2) && abs(deltaRange(i2))<= deltaLim && (abs(kRange(i4)) <= kLim  && abs(omegaRange(i3)) <=omegaLim) )
                    sinogramRebinFFTImages(i3 , i4 , i1) = sinogramRebinFFTImages(i3 , i4 , i1) + sinogramFFTImage(i3 , i4);
                    weightForSino(i3,i4,i1) = weightForSino(i3,i4,i1) + 1;
                end
            end
        end
    end
end
clear i1 i2 i3 i4 sinogramImages;

% Norm
inverseWeight = 1./weightForSino;
inverseWeight(isinf(inverseWeight)) = 0;
sinogramRebinFFTImages = sinogramRebinFFTImages.*inverseWeight;

sinogramRebinImages = zeros(length(rRange), length(phiRange), length(zDirectionRange));
for i = 1 : size(sinogramFFTImagesFinal , 3)  % sweep all z
    sinogramRebinImages(: , : , i) = ifftshift(sinogramFFTImagesFinal(: , : , i));
    sinogramRebinImages(: , : , i) = ifft2(sinogramFFTImagesFinal(: , : , i));
end
 
% Reconstruct the image
phiRange = gpuArray(phiRange);
for i =1: zDirectionNumber
    sinogramImage = gpuArray(abs(sinogramRebinImages(: , : , i )));
    sinogramImage = iradon(sinogramImage,phiRange,'Hamming',0.2,length(rRange));
    reconImage(:,:,i) = gather(sinogramImage);
end

% plot(abs(squeeze(sum(sum(reconImage,2),1))));
% plot(abs(squeeze(sum(sum(reconImage,2),3))));
% plot(abs(squeeze(sum(sum(reconImage,3),1))));
analyseXDirection=abs(squeeze(sum(sum(reconImage(790:870,:,25:37),3),1)));
plot(analyseXDirection);
[analyseXValue,analyseXIdx] = max(analyseXDirection);
analyseXInterp=spline(1:length(analyseXDirection),analyseXDirection);
analyseXFunc = @(x) ppval(analyseXInterp,x)- analyseXValue/2;
analyseXLeft = fsolve(analyseXFunc,analyseXIdx*0.99);
analyseXRight = fsolve(analyseXFunc,analyseXIdx*1.01);

analyseYDirection=abs(squeeze(sum(sum(reconImage(:,790:870,25:37),3),2)));
plot(analyseYDirection);
[analyseYValue,analyseYIdx] = max(analyseYDirection);
analyseYInterp=spline(1:length(analyseYDirection),analyseYDirection);
analyseYFunc = @(x) ppval(analyseYInterp,x)- analyseYValue/2;
analyseYLeft = fsolve(analyseYFunc,analyseYIdx*0.99);
analyseYRight = fsolve(analyseYFunc,analyseYIdx*1.01);

analyseZDirection=abs(squeeze(sum(sum(reconImage(790:870,790:870,:),1),2)));
plot(analyseZDirection);
[analyseZValue,analyseZIdx] = max(analyseZDirection);
analyseZInterp=spline(1:length(analyseZDirection),analyseZDirection);
analyseZFunc = @(x) ppval(analyseZInterp,x)- analyseZValue/2;
analyseZLeft = fsolve(analyseZFunc,analyseZIdx*0.99);
analyseZRight = fsolve(analyseZFunc,analyseZIdx*1.01);

analyseShape(:,kkk) = [(analyseXRight - analyseXLeft)*(rRange(2)-rRange(1));
                       (analyseYRight - analyseYLeft)*(rRange(2)-rRange(1));
                       (analyseZRight - analyseZLeft)*(zDirectionRange(2)-zDirectionRange(1))];
disp(['coin5.dat']);
disp(analyseShape);