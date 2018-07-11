clc; close all; clear all;
%%%%%%%%% this directory should be properly changed.
%%%%%% make sure that RANSAC and GP_toolbox had been added the path. 
% addpath('GPSegmentation\Ransac\');
% addpath('iGPR/gpml-matlab-v4.0-2016-10-19'); 
% run('iGPR/gpml-matlab-v4.0-2016-10-19/setpath.m'); 
% load data_straightRoad.mat
load data_intersection.mat
cloud = pcdownsample(pointCloud(data'), 'gridAverage', 0.1);
data = cloud.Location';
Dist = sqrt(sum(data(1:2, :).^2));
data = data(:, Dist >= 3.0 & Dist <= 60.0);
GapThr = [0.2 inf];
IS_SHOW = 1;
tic
[GrdIdx, ObsIdx, UnkownIdx, GrdPts] = GPSegFun(data, GapThr, IS_SHOW);
toc
%% Run DBSCAN Clustering Algorithm
obsData = data(:, ObsIdx); 
X = obsData(1:2, :)';
epsilon = 1.0; % the unit is meter. 
MinPts = 10;
tic
IDX = DBSCAN(X,epsilon,MinPts);
toc
%% Plot Results
figure;
hold on;
grid on;
axis equal;
xlabel('X/m');
ylabel('Y/m');
% PlotClusterinResult(X, IDX);
Colors = hsv(20);
nDim = 3; 
Legends = {};
%%%%% IDX == 0 denotes noise. 
K=max(IDX);
for i=1:K
    Tmp = obsData(:, IDX==i); 
    if length(Tmp) <= 20 %%%% discarding small obstacle
        continue; 
    end
    minRange = min(Tmp'); 
    maxRange = max(Tmp'); 
    xRange = [minRange(1) maxRange(1)]; 
    yRange = [minRange(2) maxRange(2)]; 
    zRange = [minRange(3) maxRange(3)]; 
    boxPts = [xRange([1 1 2 2 1]); 
              yRange([1 2 2 1 1]); 
              zRange([1 1 1 1 1])]; 
    id = mod(i, length(Colors)); 
    if id == 0
        id = length(Colors); 
    end
    Color = Colors(id,:);
    plot3(boxPts(1, :), boxPts(2, :), boxPts(3, :), 'color', Color, 'linestyle', '-', ...
        'linewidth', 2, 'marker', 's'); 
    if ~isempty(Tmp)
        str = sprintf('Id=%02d, Len=%02d', i-1, length(Tmp));
        pt = Tmp(:, 1);
        text(pt(1), pt(2), pt(3), str, 'FontSize', 12, 'color', 'k' );
        pcshow(Tmp', Color);
    end
end
title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);
