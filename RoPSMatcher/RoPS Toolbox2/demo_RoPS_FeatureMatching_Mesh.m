%  demo_RoPS_FeatureMatching_Mesh.m
%  Author: Yulan Guo {yulan.guo@nudt.edu.cn}
%  NUDT, China & CSSE, UWA, Australia
% This function performs feature matching on two input meshes to obtain feature
% correspondences

close all;
clc;
clear all;

keypntNum = 1000;

load ../../StanfordData/Bunny
SelIdx = [ 1 3];
pcData0 = bunny{SelIdx(1)}';
tmp = pcdownsample( pointCloud(pcData0'), 'gridAverage', 0.0015 );
pcData0 = tmp.Location';
pcData1 = bunny{SelIdx(2)}';
tmp = pcdownsample( pointCloud(pcData1'), 'gridAverage', 0.0015 );
pcData1 = tmp.Location';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load data\pointcloud_view1;
pointcloud = pcData0';
% %============================transform a pointcloud into a triangular mesh============================%
mesh = pointCloud2mesh(pointcloud,[0 0 1],0.4);                                         %other methods can also be used to perform triangulation
%============================preprocessing============================%
out = preprocessingFunc(mesh);
mesh.faceCenter = out.centroid;
mesh.faceArea = out.area;
mesh.res = out.res ;
%============================detect keypoints============================%
%keypoints are randomly seleted in this demo, any other 3D keypoint detection methods can be used
temp = randperm(length(mesh.vertices));
mesh.keypntIdx = temp(1:keypntNum);
%============================extract RoPS features at the keypoints on a mesh============================%
para.RoPS_nbSize = 15*mesh.res;
para.RoPS_binSize = 5;
para.RoPS_rotaSize = 3;
mesh.LRF =  LRFforMeshFunc(mesh, mesh.keypntIdx, para.RoPS_nbSize);
disp('LRFs calculated');  
RoPS = RoPSFunc(mesh, para.RoPS_nbSize, para.RoPS_binSize, para.RoPS_rotaSize,mesh.LRF);
mesh.RoPS = RoPS;
disp(['RoPS features generated']);  
mesh1 = mesh;
mesh1Features = [];
for keypntIdx = 1:keypntNum
    temp = trans2Dto1DFunc(mesh.RoPS{keypntIdx});
    mesh1Features = [mesh1Features; temp];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load data\pointcloud_view2;
pointcloud = pcData1';
% %============================transform a pointcloud into a triangular mesh============================%
mesh = pointCloud2mesh(pointcloud,[0 0 1],0.4);                                         %other methods can also be used to perform triangulation
%============================preprocessing============================%
out = preprocessingFunc(mesh);
mesh.faceCenter = out.centroid;
mesh.faceArea = out.area;
mesh.res = out.res ;
%============================detect keypoints============================%
%keypoints are randomly seleted in this demo, any other 3D keypoint detection methods can be used
temp = randperm(length(mesh.vertices));
mesh.keypntIdx = temp(1:keypntNum);
%============================extract RoPS features at the keypoints on a mesh============================%
para.RoPS_nbSize = 15*mesh.res;
para.RoPS_binSize = 5;
para.RoPS_rotaSize = 3;
mesh.LRF =  LRFforMeshFunc(mesh, mesh.keypntIdx, para.RoPS_nbSize);
disp('LRFs calculated');  
RoPS = RoPSFunc(mesh, para.RoPS_nbSize, para.RoPS_binSize, para.RoPS_rotaSize,mesh.LRF);
mesh.RoPS = RoPS;
disp(['RoPS features generated']);  
mesh2 = mesh;
mesh2Features = [];
for keypntIdx = 1:keypntNum
    temp = trans2Dto1DFunc(mesh.RoPS{keypntIdx});
    mesh2Features = [mesh2Features; temp];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %============================feature matching============================%
NNDRthreshold = 0.9;
corNum = 0;
kdtreeMesh1Features = KDTreeSearcher(mesh1Features,'Distance','euclidean');
for keypntIdx1 = 1:size(mesh2Features,1)
    [idxSort,distSort] = knnsearch(kdtreeMesh1Features, mesh2Features(keypntIdx1,:),'k',2,'Distance','euclidean');
    IDX = idxSort(1);
    if distSort(1)/distSort(2)<=NNDRthreshold
        corNum = corNum+1;
        corPntIdx(corNum,:) = [IDX, keypntIdx1];
        featureDis(corNum) = distSort(1);
    end
end
showCorresFunc(mesh1, mesh2, mesh1.keypntIdx(corPntIdx(:,1)), mesh2.keypntIdx(corPntIdx(:,2)), [0,200,0]);

%============================links============================%
%we may find more test datasets via the following links
% url = 'https://sites.google.com/site/yulanguo66/research-resources/3d-object-recognition-datasets';
% web(url,'-browser')


