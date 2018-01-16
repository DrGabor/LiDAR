function RoPSs = RoPSFunc(mesh,neighborSize,binSize, rotaSize, LocalFrames)

%  Author: Yulan Guo {yulan.guo@nudt.edu.cn}
%  NUDT, China & CSSE, UWA, Australia
%
% This function takes a mesh and local reference frames (LRFs) for a set of
% keypoints as an input, and generate RoPS feature descriptors for the
% keypoints as an output. An illustration is shown in Fig4 in Ref.[1].
%
% Arguments : mesh - with vertices and faces           
%                       neighborSize - the size of neighborhood to define a
%                                                  local surface for a selected keypoint
%                       binSize - number of bins for 2D plane partition
%                       rotaSize - number of rotations around each coordinate axis
%                       LocalFrames - local reference frames corresponding to all
%                                                keypoints of the input mesh

% Return :         RoPS - RoPS feature descriptors corresponding to all
%                                     keypoints
%
% Copyright : This code is written by Yulan Guo {yulan.guo@nudt.edu.cn}, NUDT. 
%               The code may be used, modified and distributed for research purposes with
%              acknowledgement of the author and inclusion this copyright information.
%References:
%           [1] Yulan Guo, Ferdous Sohel, Mohammed Bennamoun, Min Lu, Jianwei Wan. 
%           Rotational Projection Statistics for 3D Local Surface Description and Object Recognition. 
%           Internation Journal of Computer Vision. 2013, 105 (1), 63-86
%           [2] Yulan Guo, Mohammed Bennamoun, Ferdous A Sohel, Min Lu, Jianwei Wan. 
%           3D Object Recognition in Cluttered Scenes with Local Surface Features: A Survey. 
%           IEEE Transactions on Pattern Analysis and Machine Intelligence,2014
%
% Disclaimer: This code is provided as is without any warranty.

%%%%%%%%%%%% parallel computing to accelerate.
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)   % if exists no parallel computing handle, create one.
    poolsize = parpool;
    % disp('Create a Parallel Computing Pool to accelarate process' );
else
    poolsize = p.NumWorkers;
    % disp(sprintf( 'Parpool has %d workers', poolsize) );
end

keypntIdx = mesh.keypntIdx;
% BucketSize = floor(length(mesh.vertices)/100);
% kdtreeVertices = KDTreeSearcher(mesh.vertices,'Distance','euclidean','BucketSize',BucketSize);
% BucketSize = floor(length(mesh.vertices)/100);
kdtreeVertices = KDTreeSearcher(mesh.vertices,'Distance','euclidean');

interval = pi/2/rotaSize;
parfor i = 1:length(keypntIdx)
    %obtain the neighboring point of a keypoint
    keypnt = mesh.vertices(keypntIdx(i),:);
    [neighborIdx2,neighborDis]= rangesearch(kdtreeVertices,keypnt,neighborSize);    
    neighborIdx2 = cell2mat(neighborIdx2);
    neighborIdx = neighborIdx2(2:end);
    neighbNum = length(neighborIdx);
    if neighbNum<=1
        RoPSs{i,1} = ones(rotaSize*45,1)/sum(ones(rotaSize*45,1));
        continue;
    end

    %transfom the neighboring point to the local reference frame (LRF) 
    rotation = LocalFrames{i};
    neighbor = [];
    for j=1:neighbNum
        neighbor(j,:) = (mesh.vertices(neighborIdx(j),:)-mesh.vertices(keypntIdx(i),:))*inv(rotation);
    end
    
    RoPS = [];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate the sub-feature of the keypoint along the Z axis
    for rotaIdx = 1:rotaSize
        rotaAngle = (rotaIdx-1)*interval + interval/2;
        R = [cos(rotaAngle) sin(rotaAngle) 0; -sin(rotaAngle) cos(rotaAngle) 0; 0 0 1]';
        rotaNeighbor = neighbor*R;
        %projection on the XY plane 
        projNeighborXY = [rotaNeighbor(:,1),rotaNeighbor(:,2)];
        histTemp = subRoPSFunc(projNeighborXY,binSize);
        RoPS = [RoPS,histTemp];       
        %projection on the XZ plane 
        projNeighborXZ = [rotaNeighbor(:,1),rotaNeighbor(:,3)];
        histTemp = subRoPSFunc(projNeighborXZ,binSize);
        RoPS = [RoPS,histTemp];
        %projection on the YZ plane 
        projNeighborYZ = [rotaNeighbor(:,2),rotaNeighbor(:,3)];
        histTemp = subRoPSFunc(projNeighborYZ,binSize);
        RoPS = [RoPS,histTemp];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate the sub-feature of the keypoint along the Y axis
    for rotaIdx = 1:rotaSize
        rotaAngle = (rotaIdx-1)*interval + interval/2;
        R = [cos(rotaAngle) 0 sin(rotaAngle); 0 1 0 ;-sin(rotaAngle) 0 cos(rotaAngle)]'; 
        rotaNeighbor = neighbor*R;      
        %projection on the XY plane 
        projNeighborXY = [rotaNeighbor(:,1),rotaNeighbor(:,2)];
        histTemp = subRoPSFunc(projNeighborXY,binSize);
        RoPS = [RoPS,histTemp];       
        %projection on the XZ plane 
        projNeighborXZ = [rotaNeighbor(:,1),rotaNeighbor(:,3)];
        histTemp = subRoPSFunc(projNeighborXZ,binSize);
        RoPS = [RoPS,histTemp];
        %projection on the YZ plane 
        projNeighborYZ = [rotaNeighbor(:,2),rotaNeighbor(:,3)];
        histTemp = subRoPSFunc(projNeighborYZ,binSize);
        RoPS = [RoPS,histTemp];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate the sub-feature of the keypoint along the X axis
    for rotaIdx = 1:rotaSize
        rotaAngle = (rotaIdx-1)*interval + interval/2;
        R = [1,0,0; 0,cos(rotaAngle),sin(rotaAngle); 0, -sin(rotaAngle), cos(rotaAngle)]';
        rotaNeighbor = neighbor*R;
        %projection on the XY plane 
        projNeighborXY = [rotaNeighbor(:,1),rotaNeighbor(:,2)];
        histTemp = subRoPSFunc(projNeighborXY,binSize);
        RoPS = [RoPS,histTemp];       
        %projection on the XZ plane 
        projNeighborXZ = [rotaNeighbor(:,1),rotaNeighbor(:,3)];
        histTemp = subRoPSFunc(projNeighborXZ,binSize);
        RoPS = [RoPS,histTemp];
        %projection on the YZ plane 
        projNeighborYZ = [rotaNeighbor(:,2),rotaNeighbor(:,3)];
        histTemp = subRoPSFunc(projNeighborYZ,binSize);
        RoPS = [RoPS,histTemp];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RoPSs{i,1} = (RoPS)'/sum(RoPS);
end





