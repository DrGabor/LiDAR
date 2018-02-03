clc; close all; clear all;
%%%%%%%%% this directory should be properly changed. 
addpath('GPSegmentation\Ransac\');  
load data.mat
GapThr = [0.1 5.0]; 
IS_SHOW = 1; 
[EffIdx, NffIdx, UnkownIdx] = GPSegFun(data, GapThr, IS_SHOW); 