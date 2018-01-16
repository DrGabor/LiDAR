clc; close all; clear all; 
DataDir = 'StarFarming.jpg'; 
A = imread(DataDir); 
figure; 
imshow(A, []); 
tic
C = IDT(A); 
toc
figure; 
imshow(C, []); 