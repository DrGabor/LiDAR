function [ pointcloud ] = RPoS(pcData, KeyIdx )
if nargin == 0 
    clc; close all;
    load ../../StanfordData/Bunny
    pcData = bunny{1}';
    tmp = pcdownsample( pointCloud(pcData'), 'gridAverage', 0.0015 );
    pcData = tmp.Location';
    KeyIdx = 1 : 5 : size(pcData, 2);
end
if nargin == 1
    KeyIdx = 1:1:size(pcData, 2);
end
%============================detect keypoints============================%
%keypoints are randomly seleted in this demo, any other 3D keypoint detection methods can be used
pointcloud.vertices = pcData';
keypntNum = length(KeyIdx);
pointcloud.keypntIdx = KeyIdx;

%============================preprocessing============================%
kdtreeVertices = KDTreeSearcher(pointcloud.vertices,'Distance','euclidean');
[idx,dist] = knnsearch(kdtreeVertices, pointcloud.vertices,'k',2,'Distance','euclidean');
pointcloud.res = mean(dist(:,2));

%============================show the pointcloud and its keypoints============================%
pointcloudX = pointcloud;
angle = 0;%-90;
R = [1,0,0; 0,cos(angle*pi/180),sin(angle*pi/180); 0, -sin(angle*pi/180), cos(angle*pi/180)]';  %for illustration
pointcloudX.vertices = pointcloud.vertices*R;
IS_SHOW = 0;
if IS_SHOW
    figure; plot3(pointcloudX.vertices(:,1),pointcloudX.vertices(:,2),pointcloudX.vertices(:,3), 'b.'); hold on;
    plot3(pointcloudX.vertices(pointcloudX.keypntIdx,1),pointcloudX.vertices(pointcloudX.keypntIdx,2),pointcloudX.vertices(pointcloudX.keypntIdx,3),'r.', 'MarkerSize', 20);  axis equal;
end

%============================extract RoPS features at the keypoints on a point-cloud============================%
para.RoPS_nbSize = 15*pointcloud.res;
para.RoPS_binSize = 5;
para.RoPS_rotaSize = 3;
pointcloud.LRF =  LRFforPntCldFunc(pointcloud, pointcloud.keypntIdx, para.RoPS_nbSize);
% disp('LRF calculation finished');  
RoPS = RoPSFunc(pointcloud, para.RoPS_nbSize, para.RoPS_binSize, para.RoPS_rotaSize,pointcloud.LRF);
pointcloud.RoPS = RoPS;
% disp(['RoPS feature generated']);  
bTest = 1; 
end

