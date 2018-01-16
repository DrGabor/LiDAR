function [NormalInfo] = CalNormalsFun(pcData, NeighborNum, ViewPt, IS_SHOW )
if nargin == 0
    clc; close all;
    DataDir =  'E:\Dataset\Point Cloud\ply\Bunny\reconstruction\bun_zipper.ply';
    % DataDir = 'E:\Dataset\Point Cloud\ply\Kitchen_Rf6.ply'; 
    % DataDir = 'E:\Dataset\Point Cloud\ply\Torsello\angel\Angel.ply'; 
    A = pcread(DataDir);
    A = pcdownsample(A, 'random', 0.1 );
    pcData = A.Location';
    rng(100); 
    quat = rand(1, 4); 
    R = quat2rotm(quat)'; 
    T = rand(3, 1); 
    pcData = Loc2Glo( pcData, R, T ); 
    NeighborNum = 50;
    pcData = pcData(1:2, :); 
    ViewPt = zeros(2, 1);
end
if (nargin - 1) == 0
    NeighborNum = 20;
    ViewPt = zeros(3, 1);
end
if nargin < 4
    IS_SHOW = 0; 
end
%%%%%%%%%%%% calculate normals.
Md = createns(pcData');
[NNIdx, DD] = knnsearch(Md, pcData', 'k', NeighborNum );
s = struct('normal', [], 'cov', [], 'invCov', [], 'U', [], 'eig', [], 'curve', [] ); 
Dim = length(ViewPt); 
N = size(pcData, 2); 
NormalInfo = repmat(s, 1, N ); 
epsilon = 1e-4; 
USE_SVD = 0;  
parfor i = 1 : 1 : N
    idx = NNIdx(i, :); 
    Cov = cov( pcData(:, idx)' );
    if USE_SVD
        [eigenvec eigenval] = svd(Cov);
        eigenvec = eigenvec(:, 3:-1:1); 
        tmp = diag(eigenval); 
        eigenval = diag( tmp(end:-1:1) ); 
    else
        [eigenvec eigenval] = eig(Cov);
    end
    tmp = diag(eigenval); 
    NormalInfo(i).eig = tmp;
    curve = tmp(1) / sum(tmp); 
    NormalInfo(i).curve = curve; 
    U = eigenvec; 
    NormalInfo(i).U = eigenvec;
    
    % PlotElipse3D( pcData(:, idx)' ); 
    data = pcData(:, idx); 
    % S = diag([epsilon 1.0 1.0]);
    S = eye(Dim);
    S(1) = epsilon; 
    NormalInfo(i).cov = U * S * U';
    NormalInfo(i).invCov = inv(NormalInfo(i).cov); 
%     if curve > 1e-4
%         NormalInfo(i).cov = Cov; 
%     else
%         %%%%%%%%%% reconstruct the covarience, to make the cov as a plate shape. 
%         % the scale here is important, for reconstructed cov must be similar with original cov. 
%         S = eigenval(3,3) * diag([epsilon 1.0 1.0]);  
%         NormalInfo(i).cov = U * S * U'; 
%     end
    % norm corresponds the smallest eig-vector.  
    normal = eigenvec(:, 1);
    if normal' * ( ViewPt - pcData(:, i) ) < 0
        normal = normal * -1.0;
    end
    NormalInfo(i).normal = normal;
end
if IS_SHOW
    Normals = cat(2, NormalInfo(:).normal ); 
    Colors = floor( abs(Normals) * 255.0 );
    NewCloud = pointCloud( pcData', 'Normal', Normals', 'Color', uint8(Colors') );
    figure;
    view(3);
    hold on;
    grid on;
    axis equal;
    showPointCloud(NewCloud);
    [~, DD] = knnsearch( pcData', pcData', 'k', 2 ); 
    Res = mean(DD(:, 2) ); 
    EndPts = pcData + 3.0 * Res * Normals;
    tic
    for i = 1 : 10 : size(pcData, 2)
        pt0 = pcData(:, i);
        pt1 = EndPts(:, i);
        hold on; plot3( [pt0(1) pt1(1)], [pt0(2) pt1(2)], [pt0(3) pt1(3)], 'Color', abs( Normals(:, i) ) );
    end
    toc
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Normals of Point Cloud');
    
    
    cmap = floor( 255.0 * colormap(hsv(255)) );
    Curvatures = cat(1, NormalInfo(:).curve); 
    minCurve = min(Curvatures);
    maxCurve = max(Curvatures);
    CurveGap = (maxCurve - minCurve) / size(cmap, 1);
    [~, ~, id] = histcounts( Curvatures, minCurve : CurveGap : maxCurve );
    Color = cmap(id, : );
    CurveCloud = pointCloud(pcData', 'Color', uint8(Color) );
    figure;
    colormap('hsv');
    showPointCloud(CurveCloud);
    colorbar;
    title('Curvature'); 
    figure;
    showPointCloud( pcData', 'g' );
    tmp = cat(2, NormalInfo(:).eig); 
    t = abs(tmp(3, :) ./ tmp(2, :) );
    Idx = find( t >= 2.0 & t <= 4.0 );
    hold on; plot3( pcData(1, Idx), pcData(2, Idx), pcData(3, Idx), 'mp', 'MarkerSize', 5 );
    title('Curvature-based saliency point'); 
    bTest = 1;
end
end

