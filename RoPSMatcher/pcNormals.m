function [Normals, Md, varargout] = pcNormals(pcData, varargin )
if nargin == 0
    clc; close all;
    %     load ../../StanfordData/Bunny
    %     pcData = bunny{1}';
    %     tmp = pcdownsample( pointCloud(pcData'), 'gridAverage', 0.001 );
    %     pcData = tmp.Location';
    
     DataDir = 'E:\Dataset\Point Cloud\ply\Kitchen_Rf6.ply'; 
    % DataDir = 'E:\Dataset\Point Cloud\ply\Torsello\angel\Angel.ply'; 
    A = pcread( DataDir );
    A = pcdownsample(A, 'random', 0.1 );
    pcData = A.Location';
    %     x = -2.0 : 0.05 : 2.0;
    %     [x, y] = meshgrid(x, x);
    %     pcData = [ x(:)'; y(:)'; 0.000001 * x(:)' + 0.000001 * y(:)' ];
    NeighborNum = 20;
    ViewPt = zeros(3, 1);
end
if (nargin - 1) == 0
    NeighborNum = 10;
    ViewPt = zeros(3, 1);
end
if (nargin - 1) == 2
    NeighborNum = varargin{1};
    ViewPt = varargin{2};
end
%%%%%%%%%%%% parallel computing to accelerate.
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)   % if exists no parallel computing handle, create one.
    poolsize = parpool;
    % disp('Create a Parallel Computing Pool to accelarate process' );
else
    poolsize = p.NumWorkers;
    % disp(sprintf( 'Parpool has %d workers', poolsize) );
end
%%%%%%%%%%%% calculate normals.
Md = createns(pcData');
[NNIdx, DD] = knnsearch(Md, pcData', 'k', NeighborNum );
Normals = zeros(3, size(pcData, 2) );
Curvatures = zeros( 1, size(pcData, 2) );
EigArrays = zeros( 3, size(pcData, 2) );
parfor i = 1 : 1 : size(pcData, 2)
    C = cov( pcData(:, NNIdx(i, :))' );
    [eigVec eigVal] = eig(C);
    % norm corresponds the smallest eig-vector.
    [eigArray, Idx] = sort( diag(eigVal) );
    tmp = eigVec(:, Idx(1));
    if tmp' * ( ViewPt - pcData(:, i) ) < 0
        tmp = tmp * -1.0;
    end
    Normals(:, i) = tmp;
    Curvatures(i) = eigArray(1) / ( sum(eigArray) );
    EigArrays(:, i) = eigArray;
end
Miu = mean(Curvatures);
Sigma = std(Curvatures);
Idx = find( abs(Curvatures - Miu) >= 2.0 * Sigma );
if (nargout - 2) == 2
    varargout{1} = Curvatures;
    varargout{2} = EigArrays;
end
IS_SHOW = 0;
if IS_SHOW
    Colors = floor( abs(Normals) * 255.0 );
    NewCloud = pointCloud( pcData', 'Normal', Normals', 'Color', uint8(Colors') );
    figure;
    view(3);
    hold on;
    grid on;
    axis equal;
    showPointCloud(NewCloud);
    EndPts = pcData + 0.1 * Normals;
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
    figure;
    showPointCloud( NewCloud.Location, 'g' );
    t = abs(EigArrays(3, :)./EigArrays(2, :));
    Idx = find( t >= 2.0 & t <= 4.0 );
    hold on; plot3( pcData(1, Idx), pcData(2, Idx), pcData(3, Idx), 'mp', 'MarkerSize', 5 );
    bTest = 1;
end
end

