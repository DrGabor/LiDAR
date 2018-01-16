function LFs = LRFforPntCldFunc(mesh, keypntIdx, neighborSize)

%  Author: Yulan Guo {yulan.guo@nudt.edu.cn}
%  NUDT, China & CSSE, UWA, Australia
%
% This function takes a point cloud as an input, and generate local reference frames (LRFs) for a set of
% keypoints as an output.
%
% Arguments : mesh - with vertices and faces           
%                       keypntIdx - the indices of keypoints on a mesh
%                       neighborSize - the size of neighborhood to define a
%                                                  local surface for a selected keypoint
% Return :         LRFs - local reference frames (LRFs) corresponding to all
%                                     keypoints
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
% Disclaimer : This code is provided as is without any warranty.

%%%%%%%%%%%% parallel computing to accelerate.
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)   % if exists no parallel computing handle, create one.
    poolsize = parpool;
    % disp('Create a Parallel Computing Pool to accelarate process' );
else
    poolsize = p.NumWorkers;
    % disp(sprintf( 'Parpool has %d workers', poolsize) );
end

%%%%%%
% BucketSize = floor(length(mesh.vertices)/100);
% kdtreeVertices = KDTreeSearcher(mesh.vertices,'Distance','euclidean','BucketSize',BucketSize);
kdtreeVertices = KDTreeSearcher(mesh.vertices,'Distance','euclidean');
%%%%%%
parfor i = 1:length(keypntIdx)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    keypnt = mesh.vertices(keypntIdx(i),:);
    [neighborIdx,neighborDis]= rangesearch(kdtreeVertices,keypnt,neighborSize);    
    neighborIdx = cell2mat(neighborIdx);
    neighborIdx = neighborIdx(2:end);
    neighborDis = cell2mat(neighborDis);
    neighborDis = neighborDis(2:end);
    %%%%%
    M = zeros(3,3);
    dis = 0;
    for j = 1:length(neighborIdx)
        M = M+(mesh.vertices(neighborIdx(j),:)-mesh.vertices(keypntIdx(i),:))'*(mesh.vertices(neighborIdx(j),:)-mesh.vertices(keypntIdx(i),:))*(neighborSize-neighborDis(j));
        dis = dis+(neighborSize-neighborDis(j));
    end
    M = M/dis;
    if isnan(M(1,1)) ==1  
        LFs{i,1} = eye(3,3);    
        continue;     
    end
    [V,D] = eig(M);
    lamda = [D(1,1),D(2,2),D(3,3)];
    [temp, idxMinLam] = min(lamda);
    [temp, idxMaxLam] = max(lamda);
    xtemp = V(:,idxMaxLam);
    ztemp = V(:,idxMinLam);
    xPlus = 0;
    xMinus = 0;
    zPlus = 0;
    zMinus = 0;
    for j = 1:length(neighborIdx)
        if (mesh.vertices(neighborIdx(j),:)-mesh.vertices(keypntIdx(i),:))*xtemp>0
            xPlus = xPlus+1;
        else
            xMinus = xMinus+1;
        end
        if (mesh.vertices(neighborIdx(j),:)-mesh.vertices(keypntIdx(i),:))*ztemp>0
            zPlus = zPlus+1;
        else
            zMinus = zMinus+1;
        end
    end
    if xPlus>xMinus
        xAxis = xtemp;
    else
        xAxis = -xtemp;
    end
   if zPlus>zMinus
        zAxis = ztemp;
    else
        zAxis = -ztemp;        
    end 
    yAxis = cross(zAxis,xAxis);    
    rotation = [xAxis';yAxis';zAxis'];
    LFs{i,1} = rotation;
end