% % PET reconstruct from list detector responses using Fourier rebin(FORE).
% % Written by Shiwei Zhou 06/20/2020
% % Reference: Basics of PET Imaging Physics, Chemistry, and Regulations Second Edition
% %            Page 43-46, 71-73
% % Reference function: https://www.mathworks.com/help/images/ref/iradon.html
% 
% % Open bainary files, which have three columns, representing x, y and z
% fid=fopen('coin65.dat','rb');
% detectorPositions=fread(fid,'float32');
% fclose(fid);
% detectorPositions=reshape(detectorPositions,[],3);
% % figure(1);
% % clf;
% % scatter(f1(:,1),f1(:,2));
% 
% % Change the dector positions to Line of Response.
% n = 1;
% detectorLines = zeros( size(detectorPositions,1) /2 , 4);
% for idx = 1 : 2 : size(detectorPositions,1)
%     x1 = detectorPositions(idx , 1);
%     x2 = detectorPositions(idx + 1 , 1);
%     y1 = detectorPositions(idx , 2);
%     y2 = detectorPositions(idx + 1 , 2);
%     z1 = detectorPositions(idx , 3);
%     z2 = detectorPositions(idx + 1 , 3);
%     if (x1 == x2)
%         r = x1;
%         phi = -90;
%     else
%         k = (y2 - y1) / (x2 - x1);
%         phi = atan(k) / pi * 180;  % Angle of the line
%         r = (-k * x1 + y1) / sqrt(k^2 + 1); % Distance to the origin point
%     end
%     detectorLines(n,1) = r;
%     detectorLines(n,2) = phi + 180 * (phi<0);
%     detectorLines(n,3) = (z1 + z2) / 2;                  % middle point
%     detectorLines(n,4) = (z1 - z2) / sqrt( (x1 - x2)^2 + (y1 - y2)^2 );  % tangent in axial view
%     n = n + 1;
% end
% clear n idx x1 x2 y1 y2 z1 z2 k phi r;
% % figure(1);
% % clf;
% % histogram(detectorLines(:,3),10000);
% 
% % Construct the sinogram
% rMin = -8;
% rMax = 8;
% rSize = rMax - rMin;
% rNumber = 2^10;
% rRange = linspace(rMin , rMax ,rNumber);
% omegaRange = 2 * pi * (rRange(2) - rRange(1)) * (-rNumber/2 : rNumber/2 -1) /rNumber;
% 
% phiNumber = 256;
% phiRange = linspace(0,180 - 1 / phiNumber , phiNumber);
% kRange = 2 * pi * (phiRange(2) - phiRange(1)) * (-phiNumber/2 : phiNumber/2 -1) /phiNumber;
% 
% % zDirectionCluster = unique(detectorLines(:,3));
% % zDirectionMin = min(zDirectionCluster);
% % zDirectionMax = max(zDirectionCluster);
% % zDirectionNumber = length(zDirectionCluster);
% zDirectionMin = -5;
% zDirectionMax = 5;
% zDirectionNumber = 51;
% lSize = zDirectionMax - zDirectionMin;
% zDirectionRange = linspace(zDirectionMin , zDirectionMax ,zDirectionNumber);
% zDirectionDifference = zDirectionRange(2) - zDirectionRange(1);
% 
% deltaMax = min([max(detectorLines(:,4)), max(-detectorLines(:,4)) ]);
% deltaMin = -deltaMax;
% deltaNumber = 101;
% deltaRange = linspace(deltaMin , deltaMax ,deltaNumber);
% 
% sinogramImages = zeros(length(rRange), length(phiRange), length(zDirectionRange), length(deltaRange));
% for idx = 1 : size(detectorLines,1)
%     r = detectorLines(idx , 1);
%     phi = detectorLines(idx , 2);
%     z = detectorLines(idx , 3);
%     delta = detectorLines(idx , 4);
%     if ( r > rMax || r < rMin || z > zDirectionMax || z < zDirectionMin || delta > deltaMax || delta < deltaMin)
%         continue;
%     end
%     rPosition = floor((r - rMin) / (rMax - rMin) * rNumber) +1;
%     phiPosition = floor((phi)/180*phiNumber)+ 1;
%     zPosition = floor((z - zDirectionMin) / (zDirectionMax - zDirectionMin) * zDirectionNumber) +1;
%     deltaPosition = floor((delta - deltaMin) / (deltaMax - deltaMin) * deltaNumber) +1;
%     sinogramImages(rPosition , phiPosition , zPosition , deltaPosition) = sinogramImages(rPosition , phiPosition , zPosition , deltaPosition) + 1;
% end
% clear idx r phi z delta rPosition phiPosition zPosition deltaPosition;
% 
% % Fourier rebin
% deltaLim = 1;
% kLim = 2;
% rOmega = 2;
% omegaLim = 2 * (omegaRange(2) - omegaRange(1));
% % sinogramFFTImages = fft2(sinogramImages);
% % sinogramFFTImages = fftshift(sinogramFFTImages , 1);
% % sinogramFFTImages = fftshift(sinogramFFTImages , 2);
% % [omegaRangeMatrix, kRangeMatrix, zDirectionRangeMatrix, deltaRangeMatrix] = ndgrid(omegaRange, kRange, zDirectionRange, deltaRange)
% % zShiftMatrix = zDirectionRangeMatrix - deltaRangeMatrix .* kRangeMatrix ./ omegaRangeMatrix;
% % zShiftUpMatrix = ceil((zShiftMatrix - zDirectionMin) / zDirectionDifference ) + 1;
% % zShiftDownMatrix = floor((zShiftMatrix - zDirectionMin) / zDirectionDifference ) + 1;
% sinogramFFTImagesFinal = zeros(length(rRange), length(phiRange), length(zDirectionRange));
% for i1 = 1 : size(sinogramImages , 3)  % sweep all z
%     for i2 = 1 : size(sinogramImages , 4) % sweep all delta
%         sinogramFFTImage = fftshift(fft2(sinogramImages(: , : , i1 , i2)));
%         for i3 = 1 : size(sinogramImages , 1 )   % sweep all omega
%             for i4 = 1 : size(sinogramImages , 2 )   % sweep all k
%                 if(abs(kRange(i4) / omegaRange(i3))<rOmega && abs(kRange(i4))>kLim  )
%                     zShift = zDirectionRange(i1) - deltaRange(i2) * kRange(i4) / omegaRange(i3);
%                     zShiftUpIndex = ceil((zShift - zDirectionMin) / zDirectionDifference ) + 1;
%                     zShiftDownIndex = floor((zShift - zDirectionMin) / zDirectionDifference ) + 1;
%                     sinogramFFTImagesFinal(i3 , i4 , zShiftDownIndex) = sinogramFFTImagesFinal(i3 , i4 , zShiftDownIndex) ...
%                                                                     + (zShift - zDirectionRange(zShiftDownIndex)) * sinogramFFTImage(i3 , i4);
%                     sinogramFFTImagesFinal(i3 , i4 , zShiftUpIndex) = sinogramFFTImagesFinal(i3 , i4 , zShiftUpIndex) ...
%                                                                     + (zDirectionRange(zShiftUpIndex) - zShift) * sinogramFFTImage(i3 , i4);
%                 end                
%             end
%         end
%         for i3 = 1 : size(sinogramImages , 1 )   % sweep all omega
%             for i4 = 1 : size(sinogramImages , 2 )   % sweep all k        
%                 if ( abs(deltaRange(i2))<= deltaLim && abs(kRange(i4)) <= kLim  && abs(omegaRange(i3)) <=omegaLim )
%                     sinogramFFTImagesFinal(i3 , i4 , i1) = sinogramFFTImagesFinal(i3 , i4 , i1) + sinogramFFTImage(i3 , i4);
%                 end
%             end
%         end
%     end
% end
% clear i1 i2 i3 i4;

% Norm
for i1 = 1 : size(sinogramFFTImagesFinal , 3)  % sweep all z
    for i3 = 1 : size(sinogramFFTImagesFinal , 1)  % sweep all omega
        for i4 = 1 : size(sinogramFFTImagesFinal , 2)  % sweep all k
            
            if ( deltaRange(i2)<= deltaLim)
                delta2 = min([ (lSize/2 - abs(zDirectionRange(i1)))/rSize , deltaLim]);
                if (delta2 ==0)
                    sinogramFFTImagesFinal(i3 , i4 ,: ) = 0;
                else
                    sinogramFFTImagesFinal(i3 , i4 ,: ) = sinogramFFTImagesFinal(i3 , i4 ,: ) / delta2;
                end
            elseif(abs(kRange(i4) / omegaRange(i3))<rOmega && abs(kRange(i4))>kLim )
                delta1 = min([(lSize/2 - zDirectionRange(i1)) / (rSize + omegaRange(i3)/i4) , ...
                                (lSize/2 + zDirectionRange(i1)) / (rSize - omegaRange(i3)/i4) , ...
                                deltaMax - deltaMin]);
                if (delta1 ==0)
                    sinogramFFTImagesFinal(i3 , i4 ,: ) = 0;
                else
                    sinogramFFTImagesFinal(i3 , i4 ,: ) = sinogramFFTImagesFinal(i3 , i4 ,: ) / delta1;
                end
            end
            
        end
    end
end

sinogramRebinImages = zeros(length(rRange), length(phiRange), length(zDirectionRange));
for i = 1 : size(sinogramFFTImagesFinal , 3)  % sweep all z
    sinogramRebinImages(: , : , i) = ifft2(sinogramFFTImagesFinal(: , : , i));
end
 
% Reconstruct the image
phiRange = gpuArray(phiRange);
sinogramImage = gpuArray(sinogramRebinImages(: , : , end/2));
reconImage = iradon(sinogramImage,phiRange);
figure(5);
clf;
imshow(reconImage,[]);
title('Reconstructed image');
