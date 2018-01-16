clc; close all; 
nFrm = 1002; % 1000 ~ 1020 
str = sprintf('%s%04d.txt', 'R', nFrm);
DataDir = fullfile('data', str);
dataL = load(DataDir)';
options.PtsNumThr = 40;
%%%% epsilon and MinPts is DBSCAN's parameter. 
options.epsilon = 1.0;
options.MinPts = 5;
%%%% tData and tDist is GPR's parameter. 
tData = 3.0;
tDist = 0.5;
options.tData = tData;
options.tDist = tDist;
options.IS_SHOW = 1; 
%%%%%%%% iGPR to obtain valid curb points.
[GP, EffData] = iGPRFun(dataL, options);