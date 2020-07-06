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

% Change the dector positions to Line of Response.
n = 1;
detectorLines = zeros( size(detectorPositions,1) , 4);
for idx = 1  : size(detectorPositions,1)
    x1 = detectorPositions(idx , 1);
    x2 = detectorPositions(idx , 4);
    y1 = detectorPositions(idx , 2);
    y2 = detectorPositions(idx , 5);
    z1 = detectorPositions(idx , 3);
    z2 = detectorPositions(idx , 6);
    if (x1 == x2)
        r = x1;
        phi = -90;
    else
        k = (y2 - y1) / (x2 - x1);
        phi = atan(k) / pi * 180;  % Angle of the line
        r = (-k * x1 + y1) / sqrt(k^2 + 1); % Distance to the origin point
    end
    detectorLines(n,1) = r;
    detectorLines(n,2) = phi + 180 * (phi<0);
    detectorLines(n,3) = (z1 + z2) / 2;                  % middle point
    detectorLines(n,4) = (z1 - z2) / sqrt( (x1 - x2)^2 + (y1 - y2)^2 );  % tangent in axial view
    n = n + 1;
end
clear n idx x1 x2 y1 y2 z1 z2 k phi r;
% figure(1);
% clf;
% histogram(detectorLines(:,3),10000);

% Construct the sinogram
rMin = -80.5;
rMax = 80.5;
rSize = rMax - rMin;
rNumber = 2^11;
rNumber = 2334;
rRange = linspace(rMin , rMax ,rNumber);
omegaRange = 2 * pi * 1 / (rRange(2) - rRange(1)) * (-rNumber/2 : rNumber/2 -1) /rNumber;

phiNumber = 256;
phiRange = linspace(0,180, phiNumber);
kRange = 2 * pi * 1 / (phiRange(2) - phiRange(1)) * (-phiNumber/2 : phiNumber/2 -1) /phiNumber;

% zDirectionCluster = unique(detectorLines(:,3));
% zDirectionMin = min(zDirectionCluster);
% zDirectionMax = max(zDirectionCluster);
% zDirectionNumber = length(zDirectionCluster);
zDirectionMin = -24;
zDirectionMax = 24;
zDirectionNumber = 61;
lSize = zDirectionMax - zDirectionMin;
zDirectionRange = linspace(zDirectionMin , zDirectionMax ,zDirectionNumber);
zDirectionDifference = zDirectionRange(2) - zDirectionRange(1);

deltaMax = min([max(detectorLines(:,4)), max(-detectorLines(:,4)) ]);
deltaMin = -deltaMax;
deltaNumber = 81;
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
deltaLim = 0.05;
kLim = 2;
rOmega = 20;
omegaLim = 2 * (omegaRange(2) - omegaRange(1));
% sinogramFFTImages = fft2(sinogramImages);
% sinogramFFTImages = fftshift(sinogramFFTImages , 1);
% sinogramFFTImages = fftshift(sinogramFFTImages , 2);
% [omegaRangeMatrix, kRangeMatrix, zDirectionRangeMatrix, deltaRangeMatrix] = ndgrid(omegaRange, kRange, zDirectionRange, deltaRange)
% zShiftMatrix = zDirectionRangeMatrix - deltaRangeMatrix .* kRangeMatrix ./ omegaRangeMatrix;
% zShiftUpMatrix = ceil((zShiftMatrix - zDirectionMin) / zDirectionDifference ) + 1;
% zShiftDownMatrix = floor((zShiftMatrix - zDirectionMin) / zDirectionDifference ) + 1;
sinogramFFTImagesFinal = zeros(length(rRange), length(phiRange), length(zDirectionRange));
tic
for i1 = 1 : size(sinogramImages , 3)  % sweep all z
    for i2 = 1 : size(sinogramImages , 4) % sweep all delta
        if ( sum(abs(sinogramImages(: , : , i1 , i2)),'all') ==0)       % for fast calc, skip this one
            continue;
        end
        sinogramFFTImage = fftshift(fft2(sinogramImages(: , : , i1 , i2)));
        for i3 = 1 : size(sinogramImages , 1 )   % sweep all omega
            for i4 = 1 : size(sinogramImages , 2 )   % sweep all k
                if(abs(kRange(i4) / omegaRange(i3))<rOmega && (abs(kRange(i4))>kLim ||  abs(omegaRange(i3)) >omegaLim) )
                    zShift = zDirectionRange(i1)+zDirectionDifference/2 - deltaRange(i2) * kRange(i4) / omegaRange(i3);
                    zShiftUpIndex = floor((zShift - zDirectionMin) / zDirectionDifference ) + 1;
                    zShiftDownIndex = floor((zShift - zDirectionMin) / zDirectionDifference ) ;
                    if(zShiftDownIndex<zDirectionNumber+1 &&zShiftDownIndex>0)
                        sinogramFFTImagesFinal(i3 , i4 , zShiftDownIndex) = sinogramFFTImagesFinal(i3 , i4 , zShiftDownIndex) ...
                                                                    + (zShift - zDirectionRange(zShiftDownIndex)-zDirectionDifference/2) * sinogramFFTImage(i3 , i4);
                    end
                    if(zShiftUpIndex<zDirectionNumber+1 && zShiftUpIndex>0)
                        sinogramFFTImagesFinal(i3 , i4 , zShiftUpIndex) = sinogramFFTImagesFinal(i3 , i4 , zShiftUpIndex) ...
                                                                    + (zDirectionRange(zShiftUpIndex)+zDirectionDifference/2 - zShift) * sinogramFFTImage(i3 , i4);
                    end
                end
            end
        end
        for i3 = 1 : size(sinogramImages , 1 )   % sweep all omega
            for i4 = 1 : size(sinogramImages , 2 )   % sweep all k        
                if ( abs(deltaRange(i2))<= deltaLim && abs(kRange(i4)) <= kLim  && abs(omegaRange(i3)) <=omegaLim )
                    sinogramFFTImagesFinal(i3 , i4 , i1) = sinogramFFTImagesFinal(i3 , i4 , i1) + sinogramFFTImage(i3 , i4);
                end
            end
        end
    end
end
toc
save('rebin_result.mat','sinogramFFTImagesFinal');
clear i1 i2 i3 i4 sinogramImages;

% Norm
for i1 = 1 : size(sinogramFFTImagesFinal , 3)  % sweep all z
    for i3 = 1 : size(sinogramFFTImagesFinal , 1)  % sweep all omega
        for i4 = 1 : size(sinogramFFTImagesFinal , 2)  % sweep all k  
            if(abs(kRange(i4) / omegaRange(i3))<rOmega && (abs(kRange(i4))>kLim ||  abs(omegaRange(i3)) >omegaLim)  )
                delta1 = min([(lSize/2 - zDirectionRange(i1)) / (rSize + kRange(i4)/omegaRange(i3)) , ...
                                (lSize/2 + zDirectionRange(i1)) / (rSize - kRange(i4)/omegaRange(i3)) , ...
                                lSize/2/rSize]);
                if (delta1 ==0)
                    sinogramFFTImagesFinal(i3 , i4 ,i1 ) = 0;
                else
                    sinogramFFTImagesFinal(i3 , i4 ,i1 ) = sinogramFFTImagesFinal(i3 , i4 ,i1 ) / delta1;
                end
            elseif (abs(kRange(i4)) <= kLim  && abs(omegaRange(i3)) <=omegaLim )
                delta2 = min([ (lSize/2 - abs(zDirectionRange(i1)))/rSize , deltaLim]);
                if (delta2 ==0)
                    sinogramFFTImagesFinal(i3 , i4 ,i1 ) = 0;
                else
                    sinogramFFTImagesFinal(i3 , i4 ,i1 ) = sinogramFFTImagesFinal(i3 , i4 ,i1 ) / delta2;
                end
            end
            
        end
    end
end
clear i1 i3 i4;

sinogramRebinImages = zeros(length(rRange), length(phiRange), length(zDirectionRange));
for i = 1 : size(sinogramFFTImagesFinal , 3)  % sweep all z
    sinogramRebinImages(: , : , i) = ifftshift(sinogramFFTImagesFinal(: , : , i));
    sinogramRebinImages(: , : , i) = ifft2(sinogramFFTImagesFinal(: , : , i));
end
 
% Reconstruct the image
phiRange = gpuArray(phiRange);
for i =1: zDirectionNumber 
    sinogramImage = gpuArray(abs(sinogramRebinImages(: , : , i )));  
    sinogramImage = iradon(sinogramImage,phiRange,'none');
    reconImage(:,:,i) = gather(imgaussfilt(sinogramImage,2));
end

% plot(abs(squeeze(sum(sum(reconImage,2),1))));
% plot(abs(squeeze(sum(sum(reconImage,2),3))));
% plot(abs(squeeze(sum(sum(reconImage,2),1))));