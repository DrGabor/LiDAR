clc; close all; clear all; 
addpath('C:\Program Files\MATLAB\R2017a\bin\RoPSMatcher\RoPS Toolbox2\');
load data.mat
ResArray = [0.2 1.0]; 
rops0 = CalRoPSFun(data0, ResArray);
rops1 = CalRoPSFun(data1, ResArray);
tic
[R0, T0, MaxEvalNum, RefSize] = RoPSMatchFun( rops0, rops1, 1 );
toc