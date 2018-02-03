clc; close all; clear all; 
load('data.mat'); 

[~, DD] = knnsearch(Ref0', Ref0', 'k', 2); 
fineGridRes = 2*median(DD(:, 2));
xRange = -2.0:3*fineGridRes:2.0; 
yRange = xRange; 
AngArray = -180.0:2.0:180.0; 
IS_SHOW = 1; 
[ dR, dT ] = CorrelativeMatchFun(Mov0, Ref0, fineGridRes, xRange, yRange, AngArray, IS_SHOW); 