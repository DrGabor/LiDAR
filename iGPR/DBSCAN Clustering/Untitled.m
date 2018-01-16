clc; close all; 
%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
close all;

%% Load Data
DataRoot = 'F:\RG\HighWayL\'; %'F:\RG\HighWayDeparture\'; 
idx = 5178;  % 3200 
DataDir = sprintf('%sL%04d.txt', DataRoot, idx); 
dataL = load(DataDir); 
dataL = unique(dataL, 'rows'); 
DataDir = sprintf('%sR%04d.txt', DataRoot, idx); 
dataR = load(DataDir);
dataR = unique(dataR, 'rows'); 
X = [dataL; dataR]; 

%% Run DBSCAN Clustering Algorithm

epsilon=1.0;
MinPts = 5;
[ IDX, isNoise] =DBSCAN(X,epsilon,MinPts);


%% Plot Results

PlotClusterinResult(X, IDX);
title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);
ScanData = CalScanData(IDX); 

%%%%%%%%%%% ransac fit lines. 
t = 0.05; 
feedback = 1; 
[V, P, inliers] = ransacfit2Dline(X', t, feedback);
% V(1)*x+V(2)*y+V(3) = 0
x = X(inliers, 1); 
y = (-V(1)*x-V(3))/V(2); 
figure; 
hold on; 
grid on; 
axis equal; 
plot(X(:, 1), X(:, 2), 'r.' ); 
plot(X(inliers, 1), X(inliers,2), 'bo' ); 
plot(x, y, 'g.--'); 
str = sprintf('RANSAC line fit, inliers = %02d', length(inliers)); 
title(str); 
