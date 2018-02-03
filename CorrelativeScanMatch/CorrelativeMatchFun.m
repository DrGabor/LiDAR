function [ R, T ] = CorrelativeMatchFun(MovData, RefData, fineGridRes, xRange, yRange, AngArray, IS_SHOW)
%% downsample the MovData for efficiency.
Dim = size(MovData, 1); 
cCnt = 10;
coarseGridRes = cCnt * fineGridRes;
[ReduceRefData, CoorI, GridData] = VoxelizeFun( RefData, fineGridRes );
[ReduceMovData, ~, ~] = VoxelizeFun( MovData, fineGridRes );
%% generate DT image.
fineSize = max(CoorI');
coarseSize = ceil(fineSize / cCnt );
minVal = min(RefData')-1e-6;
if Dim == 3
    [X, Y, Z] = ind2sub(fineSize, 1:1:prod(fineSize));
    TGridData = bsxfun(@plus, [X; Y; Z]*fineGridRes, minVal' - ones(Dim, 1)*fineGridRes/2);
end
if Dim == 2
    [X, Y] = ind2sub(fineSize, 1:1:prod(fineSize));
    TGridData = bsxfun(@plus, [X; Y]*fineGridRes, minVal' - ones(Dim, 1)*fineGridRes/2);
end

[NNIdx, DD] = knnsearch(TGridData', ReduceRefData', 'k', 5^Dim );
Sigma = fineGridRes;
tmp = -DD(:).^2 / 2 / Sigma^2;
tmp = exp(tmp);
fineBW = zeros(fineSize);
fineBW(NNIdx(:)) = tmp;
%%%%%%%% calculate coarse BW.
[CoarseData, CoorI, ~] = VoxelizeFun( RefData, coarseGridRes );
coarseBW = zeros(coarseSize);
for id = 1 : 1 : size(CoorI, 2)
    coor = CoorI(:, id);
    tmpCoor = {};
    for i = 1 : 1 : Dim
        id0 = (coor(i)-1)*cCnt+1;
        id1 = coor(i)*cCnt;
        if id1 > fineSize(i)
            id1 = fineSize(i);
        end
        tmpCoor{end+1} = id0:1:id1;
    end
    if Dim == 3
        [X Y Z] = meshgrid(tmpCoor{1}, tmpCoor{2}, tmpCoor{3} );
        ind = sub2ind(fineSize, X(:), Y(:), Z(:) );
        coarseBW(coor(1), coor(2), coor(3)) = max( fineBW(ind) );
    end
    if Dim == 2
        [X Y] = meshgrid(tmpCoor{1}, tmpCoor{2} );
        ind = sub2ind(fineSize, X(:), Y(:) );
        coarseBW(coor(1), coor(2)) = max( fineBW(ind) );
    end
    bTest = 1;
end

%% generate fine grid.
XFineArray = xRange;
YFineArray = yRange;
[X Y] = meshgrid(XFineArray, YFineArray);
TFArray = [X(:) Y(:)]';

%% generate coarse grid.
Len = length(XFineArray);
tmpLen = floor(Len/(cCnt+1));
tmp = (cCnt+1)*(1:1:tmpLen);
A = [tmp - cCnt; tmp];
if ~isequal(A(2, end), Len)
    A(:, end+1) = [A(2, end)+1; Len];
end
A_CX = A;
XCoarseArray = ( XFineArray(A(1, :)) + XFineArray(A(2, :)) ) / 2;

Len = length(YFineArray);
tmpLen = floor(Len/(cCnt+1));
tmp = (cCnt+1)*(1:1:tmpLen);
A = [tmp - cCnt; tmp];
if ~isequal(A(2, end), Len)
    A(:, end+1) = [A(2, end)+1; Len];
end
A_CY = A;
YCoarseArray = ( YFineArray(A(1, :)) + YFineArray(A(2, :)) ) / 2;
[X Y] = meshgrid(XCoarseArray, YCoarseArray);
TCArray = [X(:) Y(:)]';
[X Y] = meshgrid(1:1:length(XCoarseArray), 1:1:length(YCoarseArray));
TCInd = [X(:) Y(:)]';
if IS_SHOW
    figure;
    hold on;
    grid on;
    axis equal;
    plot(TCArray(1, :), TCArray(2, :), 'ro');
    plot(TFArray(1, :), TFArray(2, :), 'b.');
    title('search grid');
end
% Sigma = 0.1;
%% start CSM.
Li = 1.0;
Hbest = -inf;
OptCoor = [];
LiArray = [];
HiArray = [];
HbestArray = [-inf];
TfVal = [];
parfor AngId = 1:1:length(AngArray)
    Ang = AngArray(AngId);
    R = eul2rotm(deg2rad([Ang 0 0]));
    if Dim == 2
        R = R(1:2, 1:2);
    end
    data = Loc2Glo(ReduceMovData, R', zeros(Dim, 1) );
    fVal = [];
    fValArray = [];
    for i = 1:1:size(TCArray, 2)
        T = TCArray(:, i);
        if Dim == 3
            T = [T; 0];
        end
        AftData = bsxfun(@plus, data, T);
        [Coor, EffIdx, ind] = CvtPtsToVoxelFun(AftData, minVal, coarseSize, coarseGridRes);
        fVal = sum(coarseBW(ind));
        fValArray(end+1) = fVal;
    end
    [Li, id] = max(fValArray);
    LiArray = [LiArray Li];
    %%%%% for parfor, break cannot be used....
    %     if Li < Hbest
    %         disp('break! find the optimum!');
    %         break;
    %     end
    CoarseCoor = TCInd(:, id);
    xRange = A_CX(:, CoarseCoor(1));
    yRange = A_CY(:, CoarseCoor(2));
    x = XFineArray(xRange(1):1:xRange(2));
    y = YFineArray(yRange(1):1:yRange(2));
    [x y] = meshgrid(x, y);
    tmpArray = [x(:) y(:)]';
    fValArray = [];
    for i = 1 : 1 : size(tmpArray, 2)
        T = tmpArray(:, i);
        if Dim == 3
            T = [T; 0.0];
        end
        AftData = bsxfun(@plus, data, T);
        [Coor, EffIdx, ind] = CvtPtsToVoxelFun(AftData, minVal, fineSize, fineGridRes);
        fVal = sum(fineBW(ind));
        fValArray = [fValArray fVal];
    end
    [Hi, id] = max(fValArray);
    tmp = [tmpArray(:, id); Ang; Hi];
    TfVal = [TfVal tmp];
    %     if Hi > Hbest
    %         Hbest = Hi;
    %         OptCoor = [tmpArray(:, id); Ang ];
    %     end
    HiArray = [HiArray Hi];
    bTest = 1;
end
[~, id] = max(TfVal(end, :));
OptCoor = TfVal(1:3, id);
Ang = OptCoor(3);
R = eul2rotm(deg2rad([Ang 0 0]));
T = OptCoor(1:2);
if Dim == 3
    T = [T; 0];
end
if Dim == 2
    R = R(1:2, 1:2);
end
%% visualize.
if IS_SHOW
    figure;
    hold on;
    grid on;
    plot(LiArray, 'b.');
    plot(HiArray, 'r.');
    legend('Coarse', 'Fine');
    title('fValArray');
    
    AftData = Loc2Glo(ReduceMovData, R', T);
    figure;
    hold on;
    grid on;
    axis equal;
    if Dim == 3
        view(3);
        showPointCloud(ReduceRefData', 'g');
        showPointCloud(ReduceMovData', 'r');
        showPointCloud(AftData', 'b');
    end
    if Dim == 2
        plot(ReduceRefData(1, :), ReduceRefData(2, :), 'g.');
        plot(ReduceMovData(1, :), ReduceMovData(2, :), 'r.');
        plot(AftData(1, :), AftData(2, :), 'b.');
    end
    str = sprintf('OptCoor, x = %.3f, %.3f, %.3f deg', OptCoor(1), OptCoor(2), OptCoor(3) );
    title(str);
    bTest = 1;
end
end

