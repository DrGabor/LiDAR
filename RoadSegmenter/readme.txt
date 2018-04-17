This toolbox is mainly for Gaussian Process Regression(GPR) based segmentation for autonomous driving, especially for large scenes. The algorithm is based on [1][2]. However, some major revisions are also done. 

Currently, the time cost of RoadSegmenter is around seconds, but this is not the problem if you implemented in C/C++. I tested GPR implemented with C++ using Visual Studio 2015 + i7 CPU + 16G RAM, the time cost is 20~50 ms or so.

[1] Chen, Tongtong, et al. "3D LIDAR-based ground segmentation." Pattern Recognition IEEE, 2012:8485-8490.

[2] Chen, Tongtong, et al. "Gaussian-Process-Based Real-Time Ground Segmentation for Autonomous Land Vehicles." Journal of Intelligent & Robotic Systems 76.3-4(2014):563-582.
