function [contour] = buildContour(bw, boundary, smoothing, option, value)
% USE THIS METHOD WITH bwboundary(bw,4). IT IS IMPORTANT FOR IT TO BE 4 CONNECTED.
% This methods build a contour with subpixel accuracy. In particular, it is not confused by one
% pixel thick portions of contour, 8 shapes and so on.
% The user either picks the average spacing between contour points or the number of desired contour points.

% Usage:   [ contour ] = alphaMorph( mat, alpha, contour_src, contour_dst )
%
% Arguments:
%               bw  - Binary representing the full image on which the contour lies
%         boundary  - Array of size m x 2, from cell output of C = bwboundary(bw,4), like C{1} 
%        smoothing  - Positive integer representing the smoothing factor of the contour. 
%           option  - Either nothing (default), 'Spacing' or 'Sampling'
%            value  - For 'Spacing' it represents the average distance between contour points. 
%                     For 'Sampling' it represents the exact number of desited contour points. 
% Returns:
%           contour - The contour m x 2 specified by the arguments of buildContour

% Example 1 (default):
% I = imread('.\fish\bw_fish5.png');
% BW = imbinarize(I);
% C = bwboundaries(BW,4); %IT IS IMPORTANT FOR THIS TO BE 4 CONNECTED AND ***not*** 8 CONNECTED
% C = C{1};
% contour = buildContour(BW,C); %Default: smoothing == 10, distance = 5
% 
% Example 2 ('Spacing'):
% I = imread('.\fish\bw_fish5.png');
% BW = imbinarize(I);
% C = bwboundaries(BW,4); %IT IS IMPORTANT FOR THIS TO BE 4 CONNECTED AND ***not*** 8 CONNECTED
% C = C{1};
% smoothing = 10;
% delta = 5;
% contour = buildContour(BW,C,smoothing,'Spacing',delta); 
%
% Example 3 ('Sampling'):
% I = imread('.\fish\bw_fish5.png');
% BW = imbinarize(I);
% C = bwboundaries(BW,4); %IT IS IMPORTANT FOR THIS TO BE 4 CONNECTED AND ***not*** 8 CONNECTED
% C = C{1};
% smoothing = 10;
% sampling = 200;
% contour = buildContour(BW,C,smoothing,'Sampling',sampling); 

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


%Input checking
if nargin == 2
    smoothing = 10;
    delta = 5;
    contour = build_contour_spacing(bw, boundary, smoothing, delta);
end
if nargin == 3
    error('Need to specify value for "Spacing" or "Sampling".');
end
if nargin == 4
    if(strcmp(option,'Raw'))
        contour = build_contour_raw(bw, boundary);
    else
        error('Need to specify value for "Spacing" or "Sampling".');
    end
end
if nargin == 5
    if(strcmp(option,'Spacing'))
        %check that value is positive
        contour = build_contour_spacing(bw, boundary, smoothing, value);
    elseif(strcmp(option,'Sampling'))
        %check that number is positive
        contour = build_contour_sampling(bw, boundary, smoothing, value);
    else
        error('Wrong type of strategy.');
    end
end
end

function [contour] = build_contour_spacing(bw, boundary, smoothing, value)
%Here value gives the distance between two consecutive points
%Logical initialization
bw_perim = false(size(bw));
bw_perim(sub2ind(size(bw), boundary(:,1), boundary(:,2))) = 1;
bw_down = circshift(bw,1,1);
bw_up = circshift(bw,-1,1);
bw_right = circshift(bw,1,2);
bw_left = circshift(bw,-1,2);
bw_perim_up = bw_perim & ~bw_down;
bw_perim_down = bw_perim & ~bw_up;
bw_perim_right = bw_perim & ~bw_left;
bw_perim_left = bw_perim & ~bw_right;

%Generate subpixel points
ep = 1/16;
[R,C] = find(bw_perim_up);
pts_up = [R-6*ep,C+4*ep; R-6*ep,C-4*ep];
[R,C] = find(bw_perim_down);
pts_down = [R+6*ep,C+4*ep; R+6*ep,C-4*ep];
[R,C] = find(bw_perim_right);
pts_right = [R-4*ep,C+6*ep; R+4*ep,C+6*ep];
[R,C] = find(bw_perim_left);
pts_left = [R-4*ep,C-6*ep; R+4*ep,C-6*ep];
pts = [pts_up; pts_down; pts_right; pts_left];

%Sort the points, find next closest to the current point
sorted = zeros(size(pts));
dist = pdist2(pts, pts, 'squaredeuclidean');
counter = 1;
next = 1;
for i = 1:length(dist)
    sorted(counter,:) = pts(next,:);
    dist(:,next) = 10; %no one should visit next again
    [~, ind] = min(dist(next,:));
    dist(ind, next) = 10;
    next = ind;
    counter = counter + 1;
end

%Circular smoothing
if(smoothing == 0)
    smoothed = sorted;
else
    smoothed(:,1) = conv(sorted([end-smoothing+2:end, 1:end],1),repmat(1/smoothing,1,smoothing),'valid');
    smoothed(:,2) = conv(sorted([end-smoothing+2:end, 1:end],2),repmat(1/smoothing,1,smoothing),'valid');
end

%Resample
Rc = [smoothed(:,1); smoothed(1,1)];
Cc = [smoothed(:,2); smoothed(1,2)];
dR = diff(Rc);
dC = diff(Cc);
dS = sqrt(dR.^2+dC.^2);
dS = [0; dS];
d = cumsum(dS);  % Here is the independent variable
perim = d(end);

refinement = ceil(perim / value);
delta = value;

dSi = delta*(0:refinement)';
dSi(end) = dSi(end)-.005; % Make interp1 happy
Ri = interp1(d,Rc,dSi);
Ci = interp1(d,Cc,dSi,'linear');
Ri(end)=[]; Ci(end)=[];
resampled = [Ri, Ci];


%Make sure it goes clockwise
if(oriented_area(resampled) < 0)
    flip(resampled);
end

%Store contour
contour = resampled;
end

function [contour] = build_contour_sampling(bw, boundary, smoothing, value)
%Here value gives the number of points
%Logical initialization
bw_perim = false(size(bw));
bw_perim(sub2ind(size(bw), boundary(:,1), boundary(:,2))) = 1;
bw_down = circshift(bw,1,1);
bw_up = circshift(bw,-1,1);
bw_right = circshift(bw,1,2);
bw_left = circshift(bw,-1,2);
bw_perim_up = bw_perim & ~bw_down;
bw_perim_down = bw_perim & ~bw_up;
bw_perim_right = bw_perim & ~bw_left;
bw_perim_left = bw_perim & ~bw_right;

%Generate subpixel points
ep = 1/16;
[R,C] = find(bw_perim_up);
pts_up = [R-6*ep,C+4*ep; R-6*ep,C-4*ep];
[R,C] = find(bw_perim_down);
pts_down = [R+6*ep,C+4*ep; R+6*ep,C-4*ep];
[R,C] = find(bw_perim_right);
pts_right = [R-4*ep,C+6*ep; R+4*ep,C+6*ep];
[R,C] = find(bw_perim_left);
pts_left = [R-4*ep,C-6*ep; R+4*ep,C-6*ep];
pts = [pts_up; pts_down; pts_right; pts_left];

%Sort the points, find next closest to the current point
sorted = zeros(size(pts));
dist = pdist2(pts, pts, 'squaredeuclidean');
counter = 1;
next = 1;
for i = 1:length(dist)
    sorted(counter,:) = pts(next,:);
    dist(:,next) = 10; %no one should visit next again
    [~, ind] = min(dist(next,:));
    dist(ind, next) = 10;
    next = ind;
    counter = counter + 1;
end

%Circular smoothing
if(smoothing == 0)
    smoothed = sorted;
else
    smoothed(:,1) = conv(sorted([end-smoothing+2:end, 1:end],1),repmat(1/smoothing,1,smoothing),'valid');
    smoothed(:,2) = conv(sorted([end-smoothing+2:end, 1:end],2),repmat(1/smoothing,1,smoothing),'valid');
end

%Resample
Rc = [smoothed(:,1); smoothed(1,1)];
Cc = [smoothed(:,2); smoothed(1,2)];
dR = diff(Rc);
dC = diff(Cc);
dS = sqrt(dR.^2+dC.^2);
dS = [0; dS];
d = cumsum(dS);  % independent variable
perim = d(end);

refinement = value;
delta = perim / value;

dSi = delta*(0:refinement)';
dSi(end) = dSi(end)-.005; % appease interp1
Ri = interp1(d,Rc,dSi);
Ci = interp1(d,Cc,dSi,'linear');
Ri(end)=[]; Ci(end)=[];
resampled = [Ri, Ci];

%Make sure it goes clockwise
if(oriented_area(resampled) < 0)
    flip(resampled);
end

%Store contour
contour = resampled;
end


function [contour] = build_contour_raw(bw, boundary)
%Here value gives the number of points
%Logical initialization
bw_perim = false(size(bw));
bw_perim(sub2ind(size(bw), boundary(:,1), boundary(:,2))) = 1;
bw_down = circshift(bw,1,1);
bw_up = circshift(bw,-1,1);
bw_right = circshift(bw,1,2);
bw_left = circshift(bw,-1,2);
bw_perim_up = bw_perim & ~bw_down;
bw_perim_down = bw_perim & ~bw_up;
bw_perim_right = bw_perim & ~bw_left;
bw_perim_left = bw_perim & ~bw_right;

%Generate subpixel points
ep = 1/16;
[R,C] = find(bw_perim_up);
pts_up = [R-6*ep,C+4*ep; R-6*ep,C-4*ep];
[R,C] = find(bw_perim_down);
pts_down = [R+6*ep,C+4*ep; R+6*ep,C-4*ep];
[R,C] = find(bw_perim_right);
pts_right = [R-4*ep,C+6*ep; R+4*ep,C+6*ep];
[R,C] = find(bw_perim_left);
pts_left = [R-4*ep,C-6*ep; R+4*ep,C-6*ep];
pts = [pts_up; pts_down; pts_right; pts_left];

%Sort the points, find next closest to the current point
sorted = zeros(size(pts));
dist = pdist2(pts, pts, 'squaredeuclidean');
counter = 1;
next = 1;
for i = 1:length(dist)
    sorted(counter,:) = pts(next,:);
    dist(:,next) = 10; %no one should visit next again
    [~, ind] = min(dist(next,:));
    dist(ind, next) = 10;
    next = ind;
    counter = counter + 1;
end

%Make sure it goes clockwise
if(oriented_area(sorted) < 0)
    flip(sorted);
end

%Store contour
contour = sorted;
end

function [area] = oriented_area(C)
num=size(C,1);
x0=C(1,1);
y0=C(1,2);
area=0;
c=C(2,1)-x0;
d=C(2,2)-y0;
for i=2:num-1
    a=c;
    b=d;
    c=C(i+1,1)-x0;
    d=C(i+1,2)-y0;
    area=area+a*d-c*b;
end
area=0.5*area;
end