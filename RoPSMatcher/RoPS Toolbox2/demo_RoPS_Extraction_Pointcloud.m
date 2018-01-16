%  demo_RoPS_Extraction_Pointcloud.m
%  Author: Yulan Guo {yulan.guo@nudt.edu.cn}
%  NUDT, China & CSSE, UWA, Australia
% This function extracts RoPS local features on an input point-cloud 

close all;
clc;
clear all;

load data\pointcloud_view1;

%============================detect keypoints============================%
%keypoints are randomly seleted in this demo, any other 3D keypoint detection methods can be used
pointcloud.vertices = pointcloud;
keypntNum = 100;
temp = randperm(length(pointcloud.vertices));
pointcloud.keypntIdx = temp(1:keypntNum);

%============================preprocessing============================%
kdtreeVertices = KDTreeSearcher(pointcloud.vertices,'Distance','euclidean');
[idx,dist] = knnsearch(kdtreeVertices, pointcloud.vertices,'k',2,'Distance','euclidean');
pointcloud.res = mean(dist(:,2));

%============================show the pointcloud and its keypoints============================%
pointcloudX = pointcloud;
angle = 0;%-90;
R = [1,0,0; 0,cos(angle*pi/180),sin(angle*pi/180); 0, -sin(angle*pi/180), cos(angle*pi/180)]';  %for illustration
pointcloudX.vertices = pointcloud.vertices*R;
figure; plot3(pointcloudX.vertices(:,1),pointcloudX.vertices(:,2),pointcloudX.vertices(:,3), 'b.'); hold on;
plot3(pointcloudX.vertices(pointcloudX.keypntIdx,1),pointcloudX.vertices(pointcloudX.keypntIdx,2),pointcloudX.vertices(pointcloudX.keypntIdx,3),'r.', 'MarkerSize', 20);  axis equal;

%============================extract RoPS features at the keypoints on a point-cloud============================%
para.RoPS_nbSize = 15*pointcloud.res;
para.RoPS_binSize = 5;
para.RoPS_rotaSize = 3;
pointcloud.LRF =  LRFforPntCldFunc(pointcloud, pointcloud.keypntIdx, para.RoPS_nbSize);
disp('LRF calculation finished');  
RoPS = RoPSFunc(pointcloud, para.RoPS_nbSize, para.RoPS_binSize, para.RoPS_rotaSize,pointcloud.LRF);
pointcloud.RoPS = RoPS;
disp(['RoPS feature generated']);  

%============================links============================%
%we may find more test datasets via the following links
url = 'https://sites.google.com/site/yulanguo66/research-resources/3d-object-recognition-datasets';
web(url,'-browser')



