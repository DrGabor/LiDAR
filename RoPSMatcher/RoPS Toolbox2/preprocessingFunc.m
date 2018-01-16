function out = preprocessingFunc(mesh)

%  Author: Yulan Guo {yulan.guo@nudt.edu.cn}
%  NUDT, China & CSSE, UWA, Australia

% This function takes a mesh as an input, and calculate the face centriod, face area and  the mesh resolution.
%
% Arguments :   mesh - with vertices and faces           
% Return :         out.area - the area of each face in the input triangular mesh
%                      out.centroid - the centroid of each face in the input triangular mesh
%                      out.res - the resolution of the input triangular mesh
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

vertex = mesh.vertices;
face = mesh.faces;

[vertex,face] = check_face_vertex(vertex,face);
nface = size(face,2);
nvert = size(vertex,2);
centroid = zeros(3,nface);
normalf = crossp(vertex(:,face(2,:))-vertex(:,face(1,:)), ...
                  vertex(:,face(3,:))-vertex(:,face(1,:)) );
for i=1:nface
    f = face(:,i);
    area(i,1) = 0.5*norm(normalf(:,i));
    for j=1:3
        centroid(:,i) = centroid(:,i) + vertex(:,f(j));
    end
end
centroid = centroid./3;

edges = compute_edges((mesh.faces)');
nedge = length(edges);
for i=1:nedge
   dis(i) = norm(vertex(:,edges(1,i))-vertex(:,edges(2,i)));
end
res = mean(dis);

out.area = area;
out.centroid = centroid';
out.res = res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
z = x;
z(1,:) = x(2,:).*y(3,:) - x(3,:).*y(2,:);
z(2,:) = x(3,:).*y(1,:) - x(1,:).*y(3,:);
z(3,:) = x(1,:).*y(2,:) - x(2,:).*y(1,:);
