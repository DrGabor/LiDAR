function mesh = pointCloud2mesh(data, refNormal, stdTol)

% mesh = meshD(data, refNormal, stdTol)

% Author : Ajmal Saeed Mian {ajmal@csse.uwa.edu.au}
%           Computer Science. Univ of Western Australia
%
% This function takes data points performs triangulation on it, filters out
% incorrecp polygons and outputs a mesh data structure like the newMesh
% function.
%
% Arguments : data - Nx3 vertex coordinates [x y z] of the pointcloud
%             stdTol - (optional) tolerance for edge filtering. default is 0.6
%             
%             refNormal - (optional) 1x3 vector in the sensor direction
%                         =[0 0 1] if the sensor looking towards the -z_axis
%
% Return : mesh - mesh data structure
%                       vertices: Nx3 vertex coordinates
%                       faces: M faces using index numbers of the vertices
%                       resolution: the mean edge length of faces
%                       stdeviation: the standard deviation o edge lengths
%                       triangleNormals: Mx3 normal vectors of each triangle
%                       vertexNormals: Nx3 normal vectors of each vertex
%                       vertexNtriangles: Nx1 cell of neighboring faces 
%                                           of each vertex
%                       triangleNtriangles: Mx1 cell of nieghboring faces
%                                               of each triangle
%
% Copyright : This code is written by Ajmal Saeed Mian {ajmal@csse.uwa.edu.au}
%              Computer Science, The University of Western Australia. The code
%              may be used, modified and distributed for research purposes with
%              acknowledgement of the author and inclusion this copyright information.
%
% Disclaimer : This code is provided as is without any warrantly.

warning off MATLAB:divideByZero;
if nargin == 1
    PC = princomp(data);
    data = data*PC;
    refNormal = [0 0 1];
    refNormal = refNormal * PC;
end

if nargin < 3
    stdTol = 0.6;
end

tri = delaunay(data(:,1),data(:,2));
tri(:,4) = 0; % initialize 4th column to store maximum edge length

edgeLength = [sqrt(sum((data(tri(:,1),:) - data(tri(:,2),:)).^2,2)),...
        sqrt(sum((data(tri(:,2),:) - data(tri(:,3),:)).^2,2)),...
        sqrt(sum((data(tri(:,3),:) - data(tri(:,1),:)).^2,2))];

tri(:,4) = max(edgeLength,[],2);

resolution = mean(edgeLength(:));
stdeviation = std(edgeLength(:));
filtLimit = resolution + stdTol*stdeviation;

bigTriangles = find(tri(:,4) > filtLimit); %find index numbers of triagles with edgelength more than filtLimit
tri(bigTriangles,:) = []; % remove all faces with edgelength more than filtlimit
tri(:,4) = []; % remove the max edgeLength column

edgeLength(bigTriangles,:) = []; % remove edges belonging to faces which are removed
edgeLength = edgeLength(:); 
resolution = mean(edgeLength); % find the mean of the remaining edges
stdeviation = std(edgeLength);

mesh = [];
if nargin < 2
    data = data*PC';% multiply the data points by the inverse PC
    refNormal = refNormal * PC';
end
mesh.vertices = data;  
mesh.faces = tri;
mesh.resolution = resolution;
mesh.stdeviation = stdeviation;

noOfpolygons = size(tri,1);
noOfpoints = size(data,1);
mesh.triangleNormals = zeros(noOfpolygons,3); % innitialize a matrix to store polygon normals
mesh.vertexNormals = zeros(noOfpoints,3); % innitialize a matrix to store point normals
mesh.vertexNtriangles = cell(noOfpoints, 1); %a cell array to store neighbouring polygons for the current point
mesh.triangleNtriangles = cell(noOfpolygons, 1); % to store neighbors of current polygon

for ii = 1:noOfpolygons %find normals of all polygons
    %indices of the points from which the polygon is made
    pointIndex1 = mesh.faces(ii,1);
    pointIndex2 = mesh.faces(ii,2);
    pointIndex3 = mesh.faces(ii,3);
    
    %coordinates of the points
    point1 = mesh.vertices(pointIndex1,:);
    point2 = mesh.vertices(pointIndex2,:);
    point3 = mesh.vertices(pointIndex3,:);
    
    vector1 = point2 - point1;
    vector2 = point3 - point2;
    
    normal = cross(vector1,vector2);
    normal = normal / norm(normal);
    
    theta = acos(dot(refNormal, normal));
    if theta > pi/2
        normal = normal * (-1);
        a = mesh.faces(ii,2);
        mesh.faces(ii,2) = mesh.faces(ii,1);
        mesh.faces(ii,1) = a;
    end
    
    mesh.triangleNormals(ii,:)=normal;   
            
    %make entry of this polygon as the neighbouring polygon of the three
    %vertex points    
    mesh.vertexNtriangles(pointIndex1,1)={[mesh.vertexNtriangles{pointIndex1,1} ii]};
    mesh.vertexNtriangles(pointIndex2,1)={[mesh.vertexNtriangles{pointIndex2,1} ii]};
    mesh.vertexNtriangles(pointIndex3,1)={[mesh.vertexNtriangles{pointIndex3,1} ii]};    
end

for ii = 1:noOfpoints %find normals of all points
    polys = mesh.vertexNtriangles{ii};% get neighboring polygons to this point
    normal2 = zeros(1,3);
        
    for jj = 1 : size(polys,1)
        normal2 = normal2 + mesh.triangleNormals(polys(jj),:);
    end
    
    normal2 = normal2 / norm(normal2);
    mesh.vertexNormals(ii,:) = normal2;
end

for ii = 1 : noOfpolygons % find neighbouring polygons of all polygons
    polNeighbor = [];
    for jj = 1 : 3
        polNeighbor = [polNeighbor mesh.vertexNtriangles{mesh.faces(ii,jj)}];
    end
    polNeighbor = unique(polNeighbor);
    polNeighbor = setdiff(polNeighbor, [ii]);
    mesh.triangleNtriangles(ii,1)={[polNeighbor]};
end