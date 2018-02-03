README
================================================================================================
Version 1.0, 04-Feb-2018

The code has been tested with MATLAB 2017a/b on a PC with 64-bit windows 7/10.

================================================================================================

Use of this code is free for research purposes only.

================================================================================================
Correlative scan matching(CSM)[1] is a nearly global map/scan matching algorithm for 2D data. CSM can provide a quite good initial value for other accurate matching algorithms, like ICP or NDT. CSM can also used as loop closure detector. 
The basic theory lies in two layer branch and bound. It should be noted that the original CSM is implemented in C++, GPU boosted OpenGL. This MATLAB code is just for investigating the theory in CSM.  

[1]. Olson, E.B. Real-time correlative scan matching. in Robotics and Automation, 2009. ICRA'09. IEEE International Conference on. 2009: IEEE.
