% PET reconstruct from list detector responses.
% Written by Shiwei Zhou 06/17/2020
% Reference: Basics of PET Imaging Physics, Chemistry, and Regulations Second Edition
%            Page 43-46, 71-73
% Reference function: https://www.mathworks.com/help/images/ref/iradon.html

% % Open bainary files, which have three columns, representing x, y and z
% fid=fopen('coin65.dat','rb');
% detectorPositions=fread(fid,'float32');
% fclose(fid);
% detectorPositions=reshape(detectorPositions,[],3);
% % figure(1);
% % clf;
% % scatter(f1(:,1),f1(:,2));

% Select a slice along z direction
interestedDetectorPositions=[];
interestedDetectorMinZ = -0.5;
interestedDetectorMaxZ = 0.5;
n = 1;
for idx = 1 : 2 : size(detectorPositions,1)
    if (detectorPositions(idx,3)>interestedDetectorMinZ &&...
        detectorPositions(idx,3)<interestedDetectorMaxZ &&...
        detectorPositions(idx+1,3)>interestedDetectorMinZ &&...
        detectorPositions(idx+1,3)<interestedDetectorMaxZ )
    
        interestedDetectorPositions(n,:) = detectorPositions(idx,:);
        interestedDetectorPositions(n+1,:) = detectorPositions(idx+1,:);
        n = n + 2;
    end
end
clear n idx;
% figure(2);
% clf;
% scatter3(interestedDetectorPositions(:,1),interestedDetectorPositions(:,2),interestedDetectorPositions(:,3));

% Change the dector positions to Line of Response.
n = 1;
interestedDetectorLines = zeros( size(interestedDetectorPositions,1) /2 , 2);
for idx = 1 : 2 : size(interestedDetectorPositions,1)
    x1 = interestedDetectorPositions(idx , 1);
    x2 = interestedDetectorPositions(idx + 1 , 1);
    y1 = interestedDetectorPositions(idx , 2);
    y2 = interestedDetectorPositions(idx + 1 , 2);
    if (x1 == x2)
        r = x1;
        phi = -90;
    else
        k = (y2 - y1) / (x2 - x1);
        phi = atan(k) / pi * 180;  % Angle of the line
        r = (-k * x1 + y1) / sqrt(k^2 + 1); % Distance to the origin point
    end
    interestedDetectorLines(n,1) = r;
    interestedDetectorLines(n,2) = phi;
    n = n + 1;
end
clear n idx x1 x2 y1 y2 k phi r;

% Construct the sinogram
rRange = -8:0.01:8;
phiRange = 0:0.1:179.9;
% phiRange = 0:1:179;
sinogramImage = zeros(length(rRange),length(phiRange),'gpuArray');
for idx = 1 : 2 : size(interestedDetectorLines,1)
    r = interestedDetectorLines(idx , 1);
    phi = interestedDetectorLines(idx , 2);
    rPosition = floor((r+8)/0.01) +1;
    phiPosition = floor((phi+90)/0.1)+ 1;
    sinogramImage(rPosition , phiPosition) = sinogramImage(rPosition , phiPosition) + 1;
end
figure(3);
clf;
imshow(sinogramImage,[]);
title('sinogram');
clear idx phi r rPosition phiPosition;

% Reconstruct the image
phiRange = gpuArray(phiRange);
reconImage = iradon(sinogramImage,phiRange);
figure(5);
clf;
imshow(reconImage,[]);
title('Reconstructed image');
