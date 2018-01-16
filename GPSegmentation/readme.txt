README
================================================================================================
Version 1.0, 16-Jan-2018

The code has been tested with MATLAB 2017a/b on a PC with 64-bit windows 7/10.

================================================================================================

Use of this code is free for research purposes only.

================================================================================================
For autonomous driving using 3D LiDAR, the first task is classifying point cloud as ground and obstacles. This toolbox implements Gaussian Process Regression(GPR) based segmentation, which is quite accurate and robust for real urban traffic scenes.
The algorithm is based on following papers:
[1] Chen, Tongtong, et al. "3D LIDAR-based ground segmentation." Pattern Recognition IEEE, 2012:8485-8490.
[2] Chen, Tongtong, et al. "Gaussian-Process-Based Real-Time Ground Segmentation for Autonomous Land Vehicles." Journal of Intelligent & Robotic Systems 76.3-4(2014):563-582.
