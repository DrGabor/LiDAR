For two 3D point cloud, sometimes we need align them coarsely by hand or feature matching. Rotational projection statistics(RoPS)[1] is a robust 3D feature for point cloud just like SIFT feature for image.
The "RoPS Toolbox2" is written by Prof. Yulan Guo, including RoPS calculation and RoPS feature matching. 
I add the rotation/translation calculation in RoPSMatchFun.m.
It should be noted that though RoPS feature calculation had been implemented in PCL 1.8.0, the feature matching and rotation/translation calculation are not added in PCL 1.8.0. 

[1] Guo, Yulan, et al. "Rotational Projection Statistics for 3D Local Surface Description and Object Recognition." International Journal of Computer Vision 105.1(2013):63-86.