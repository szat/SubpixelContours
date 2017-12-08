% Copyright (c) Adrian Szatmari
% Author: Adrian Szatmari
% Date: 2017-11-30
% License: MIT, patent permitting
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

I = imread('..\data\bw_fish5.png');
BW = imbinarize(I);
C = bwboundaries(BW,4); %IT IS IMPORTANT FOR THIS TO BE 4 CONNECTED AND ***not*** 8 CONNECTED
C = C{1};

%Default: smoothing == 10, distance = 5
contour = buildContour(BW,C);
%Distance between points
dist = sqrt((contour(:,1) - circshift(contour(:,1),1)).^2 + (contour(:,2) - circshift(contour(:,2),1)).^2);
figure
plot(1:length(dist),dist,'r');
title('Distance between consecutive points, option "Default"');
%Viz
figure
imshow(I);
hold on 
plot(contour(:,2),contour(:,1),'-.xr');
plot(contour(1,2),contour(1,1),'o','LineWidth',2);
hold off

%Specify the amount of space between points
smoothing = 10;
space = 5.5;
contour = buildContour(BW,C,smoothing,'Spacing',space);
%Distance between points
dist = sqrt((contour(:,1) - circshift(contour(:,1),1)).^2 + (contour(:,2) - circshift(contour(:,2),1)).^2);
figure
plot(1:length(dist),dist,'g');
title('Distance between consecutive points, option "Spacing"');
%Viz
figure
imshow(I);
hold on 
plot(contour(:,2),contour(:,1),'-.xg');
plot(contour(1,2),contour(1,1),'o','LineWidth',2);
hold off

%Specify the number of points on the contour
smoothing = 10;
nb_points = 200;
contour = buildContour(BW,C,smoothing,'Sampling',nb_points);
%Distance between points
dist = sqrt((contour(:,1) - circshift(contour(:,1),1)).^2 + (contour(:,2) - circshift(contour(:,2),1)).^2);
figure
plot(1:length(dist),dist,'b');
title('Distance between consecutive points, option "Sampling"');
%Viz
figure
imshow(I);
hold on 
plot(contour(:,2),contour(:,1),'-.bx');
plot(contour(1,2),contour(1,1),'o','LineWidth',2);
hold off

%Raw: no smoothing, no resampling, the smoothing term does nothing
contour = buildContour(BW,C,0,'Raw');
%Distance between points
dist = sqrt((contour(:,1) - circshift(contour(:,1),1)).^2 + (contour(:,2) - circshift(contour(:,2),1)).^2);
figure
plot(1:length(dist),dist,'c');
title('Distance between consecutive points, option "Raw"');
%Viz
figure
imshow(I);
hold on 
plot(contour(:,2),contour(:,1),'-.xc');
plot(contour(1,2),contour(1,1),'o','LineWidth',2);
hold off
