function [R0, T0, maxEvalNum, RefSize] = RoPSMatchFun(RoPSInfo_ref, RoPSInfo_mov, IS_SHOW )
tic;
pcData0 = RoPSInfo_ref.vertices;
pcData1 = RoPSInfo_mov.vertices;
[NNIdx, DD] = knnsearch( RoPSInfo_ref.RoPS, RoPSInfo_mov.RoPS, 'k', 2 );
%%%%%%%%%%%%%%%% find reliable features.
tmp = DD(:, 1) ./ DD(:, 2);
Idx = find( tmp <= 0.9 );
MapIdx0 = NNIdx(:, 1);      % MapIdx0 and MapIdx1 is the index of RoPS features(far few than point cloud).
MapIdx1 = 1:1:length(RoPSInfo_mov.keyIdx);
SIdx0 = RoPSInfo_ref.keyIdx(NNIdx(:, 1));   % SIdx0 and SIdx1 is the index of point cloud.
SIdx1 = RoPSInfo_mov.keyIdx;

nFeatureTime = toc;
%%%%%%%%%%%%%%%%% compute rigid transforms. Using brutal force methods.
tic
RefData  = pcData0(:, SIdx0 );
MovData  = pcData1(:, SIdx1 );
params.DistThr = RoPSInfo_ref.CoarseRes;
params.MaxIter = Inf;
DistThr = params.DistThr;
USE_ROPS_REF = 1;
if USE_ROPS_REF
    EvalNum = [];
    Tf = zeros(3, 4, size(MovData, 2) );
    parfor i = 1 : 1 : size(MovData, 2)
        LRF0 = RoPSInfo_mov.LRF{ MapIdx1(i) };   % estimates coordinate0 to coordinate1
        LRF1 = RoPSInfo_ref.LRF{ MapIdx0(i) };
        pt0  = MovData(:, i)';
        pt1  = RefData(:, i)';
        EstR = LRF0' * LRF1;
        EstT = pt1 - pt0 * EstR;
        R = EstR';
        T = EstT';
        AftData = Loc2Glo( MovData, R', T );
        tmpDist = PointsNorm( AftData - RefData );
        EvalNum = [ EvalNum length( find( tmpDist < DistThr^2 ) )];
        Tf(:, :, i) = [R T];
    end
    [maxEvalNum Idx] = max(EvalNum);
    R0 = Tf(:, 1:end-1, Idx);
    T0 = Tf(:, end, Idx );
    tmp = bsxfun( @plus, R0 * MovData, T0 ) - RefData;
    tmpDist = tmp(1, :).^2 + tmp(2, :).^2 + tmp(3, :).^2;
    Idx = find( tmpDist < DistThr^2 );
    [R0, T0, tmpError ] = RecoverRT( MovData, RefData, Idx);
else
    [R0, T0, tmpError] = ExtensiveRT( MovData, RefData, params );
end
%%%%%%%%%%%%%%% Refine the transformations.
maxEvalNum = length(find(tmpError < DistThr^2));
RefSize = size(RefData, 2); 
nRecoverTime = toc;

if IS_SHOW
    figure;
    subplot(121);
    axis equal;
    grid on;
    hold on;
    plot3( pcData0(1, :), pcData0(2, :), pcData0(3, :), 'b.');
    plot3( pcData0(1, SIdx0), pcData0(2, SIdx0), pcData0(3, SIdx0), 'r.', 'MarkerSize', 20);
    title( sprintf( 'Key Points = %04d', length(SIdx0)) );
    subplot(122);
    axis equal;
    grid on;
    hold on;
    plot3( pcData1(1, :), pcData1(2, :), pcData1(3, :), 'b.');
    plot3( pcData1(1, SIdx1), pcData1(2, SIdx1), pcData1(3, SIdx1), 'r.', 'MarkerSize', 20);
    title( sprintf( 'Key Points = %04d', length(SIdx1)) );
    
    figure;
    axis equal;
    grid on;
    hold on;
    plot3( pcData0(1, :), pcData0(2, :), pcData0(3, :), 'c.');
    plot3( pcData0(1, SIdx0), pcData0(2, SIdx0), pcData0(3, SIdx0), 'r.', 'MarkerSize', 20);
    plot3( pcData1(1, :), pcData1(2, :), pcData1(3, :), 'b.');
    plot3( pcData1(1, SIdx1), pcData1(2, SIdx1), pcData1(3, SIdx1), 'r.', 'MarkerSize', 20);
    for i = 1 : 10 : length(SIdx0)
        pt0 = pcData0( :, SIdx0(i) );
        pt1 = pcData1( :, SIdx1(i) );
        plot3( [pt0(1) pt1(1)], [pt0(2) pt1(2)], [pt0(3) pt1(3)], 'm' );
    end
    title( 'Corresponds' );
    
    %%%%% after coarse alignment.
    figure;
    grid on;
    hold on;
    view(3);
    showPointCloud(pcData0', 'g' );
    showPointCloud(pcData1', 'r' );
    AftPts = bsxfun(@plus, R0 * pcData1, T0 );
    showPointCloud(AftPts', 'b' );
    title(sprintf( 'Coarse Registration, T = %.2f %.2f %.2f', T0(1), T0(2), T0(3) ) );
end
end



