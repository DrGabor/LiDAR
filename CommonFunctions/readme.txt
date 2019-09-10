README
================================================================================================
Version 1.0, 16-Jan-2018

This package contains the MATLAB implementation of common functions for 3D point cloud processing, especially for HDL-64E S3. 

The code has been tested with MATLAB 2017a/b on a PC with 64-bit windows 7/10.

================================================================================================

Use of this code is free for research purposes only.

================================================================================================
CoorTf: 
make sure that <R, T> is coordinate transform, i.e. if Ang = sita, R = [cos(sita) sin(sita); -sin(sita) cos(sita)], which is transpose of rigid transformation.
Loc2Glo(data, R, T) local to global transformation.
Glo2Loc(data, R, T)
Loc2Loc(data, R, T)
Point Registration:
p2pICP() point to point Iterative Closest Point(ICP)
p2plICP() point to plane ICP
TrimmedICP() Trimmed ICP

Velodyne 64E-S3 Related
CalNormalsFun() calculate normals for 3D point cloud. 
