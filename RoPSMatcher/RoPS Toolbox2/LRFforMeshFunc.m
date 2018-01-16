function LRFs = LRFforMeshFunc(mesh, keypntIdx, neighborSize)

%  Author: Yulan Guo {yulan.guo@nudt.edu.cn}
%  NUDT, China & CSSE, UWA, Australia
%
% This function takes a mesh as an input, and generate local reference frames (LRFs) for a set of
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

BucketSize = floor(length(mesh.vertices)/100);
kdtreeFaceCenter = KDTreeSearcher(mesh.faceCenter,'Distance','euclidean','BucketSize',BucketSize);

for i = 1:length(keypntIdx)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    keypnt = mesh.vertices(keypntIdx(i),:);
    [neighborFaceIdx,neighborFaceDis]= rangesearch(kdtreeFaceCenter,keypnt,neighborSize);    
    neighborFaceIdx = cell2mat(neighborFaceIdx);
    neighborFaceDis = cell2mat(neighborFaceDis);
    neighborFaceLength = length(neighborFaceIdx);
    Len1 = floor(neighborFaceLength*0.8);
    neighborIdx = [];    
    M = zeros(3,3);
    totalArea = 0;
    tiotalDisW = 0;
    for j=1:neighborFaceLength
        vertIdx = mesh.faces(neighborFaceIdx(j),:);           
        dis1(j) = norm(keypnt-mesh.vertices(vertIdx(1),:));
        dis2(j) = norm(keypnt-mesh.vertices(vertIdx(2),:));
        dis3(j) = norm(keypnt-mesh.vertices(vertIdx(3),:));    
        if dis1(j)<neighborSize && dis2(j)<neighborSize && dis3(j)<neighborSize
            neighborIdx = [neighborIdx,vertIdx];
            area(j) = mesh.faceArea(neighborFaceIdx(j));
            totalArea = totalArea+area(j);            
            disW(j) = (neighborSize-neighborFaceDis(j))^2;
            tiotalDisW = tiotalDisW+disW(j);
            center_i = zeros(3,3);
            for ii=1:3
                temp11 = mesh.vertices(vertIdx(ii),:)-keypnt;
                for jj=1:3
                    center_i = center_i+temp11'*(mesh.vertices(vertIdx(jj),:)-keypnt);
                end
                center_i = center_i+temp11'*temp11;
                displace(ii,:) = temp11;
            end
            displaceTotal{j} = displace;   
            M = M+center_i*area(j)*disW(j);
        end
    end
    
    M = M/totalArea/tiotalDisW;    
    if isnan(M(1,1)) ==1  
        LRFs{i,1} = eye(3,3);    
        continue;     
    end
        
    [V,D] = eig(M);
    lamda = [D(1,1),D(2,2),D(3,3)];
    [temp, idxMinLam] = min(lamda);
    [temp, idxMaxLam] = max(lamda);
    xtemp = V(:,idxMaxLam);
    ztemp = V(:,idxMinLam);
    xDisplace = 0;
    zDisplace = 0;

    %disambiguating the sign of x- y- and z- axis
    for j=1:neighborFaceLength
        if j>=Len1
            rf = dis1(j)<neighborSize && dis2(j)<neighborSize && dis3(j)<neighborSize;
            if rf==0
                continue;
            end
        end
        displace = displaceTotal{j};   
        for ii=1:3
            xDisplace = xDisplace + displace(ii,:)*xtemp*area(j)*disW(j);
            zDisplace = zDisplace + displace(ii,:)*ztemp*area(j)*disW(j);
        end      
    end
    if xDisplace>0
        xAxis = xtemp;
    else
        xAxis = -xtemp;
    end
   if zDisplace>0
        zAxis = ztemp;
    else
        zAxis = -ztemp;        
    end 
    yAxis = cross(zAxis,xAxis);    
    %get the local coordinates of neighboring points
    rotation = [xAxis';yAxis';zAxis'];
    LRFs{i,1} = rotation;
end