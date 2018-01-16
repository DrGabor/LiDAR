%  demo_RoPS_Extraction_Mesh.m
%  Author: Yulan Guo {yulan.guo@nudt.edu.cn}
%  NUDT, China & CSSE, UWA, Australia
% This function extracts RoPS local features on a mesh (converted from an input point-cloud) 

close all;
clc;
clear all;

load data\pointcloud_view1;

% %============================transform a pointcloud into a triangular mesh============================%
mesh = pointCloud2mesh(pointcloud,[0 0 1],0.4);                                         %other methods can also be used to perform triangulation

%============================preprocessing============================%
out = preprocessingFunc(mesh);
mesh.faceCenter = out.centroid;
mesh.faceArea = out.area;
mesh.res = out.res ;

%============================detect keypoints============================%
%keypoints are randomly seleted in this demo, any other 3D keypoint detection methods can be used
keypntNum = 100;
temp = randperm(length(mesh.vertices));
mesh.keypntIdx = temp(1:keypntNum);

%============================show the mesh and its keypoints============================%
meshX = mesh;
angle = 0;%-90;
R = [1,0,0; 0,cos(angle*pi/180),sin(angle*pi/180); 0, -sin(angle*pi/180), cos(angle*pi/180)]';
meshX.vertices = mesh.vertices*R;
figure; trisurf(meshX.faces,meshX.vertices(:,1),meshX.vertices(:,2),meshX.vertices(:,3)); axis equal;
axis image; shading interp;lighting phong; view([0,0]);  hold on; colormap([1,1,1]*0.9);camlight right;hold on;axis off;
plot3(meshX.vertices(mesh.keypntIdx,1),meshX.vertices(mesh.keypntIdx,2),meshX.vertices(mesh.keypntIdx,3),'r.');

%============================extract RoPS features at the keypoints on a mesh============================%
para.RoPS_nbSize = 15*mesh.res;
para.RoPS_binSize = 5;
para.RoPS_rotaSize = 3;
mesh.LRF =  LRFforMeshFunc(mesh, mesh.keypntIdx, para.RoPS_nbSize);
disp('LRFs calculated');  
RoPS = RoPSFunc(mesh, para.RoPS_nbSize, para.RoPS_binSize, para.RoPS_rotaSize,mesh.LRF);
mesh.RoPS = RoPS;
disp(['RoPS features generated']);  

%============================links============================%
%we may find more test datasets via the following links
url = 'https://sites.google.com/site/yulanguo66/research-resources/3d-object-recognition-datasets';
web(url,'-browser')




