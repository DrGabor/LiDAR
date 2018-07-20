function [R0, T0, tmpError] = ExtensiveRT(MovData, RefData, params )
if nargin == 0
    clc; close all; 
    load ../../StanfordData/Bunny
    SelIdx = [2 3];
    CoarseRes = 0.01;
    pcData0 = bunny{SelIdx(1)}';
    R = eul2rotm( deg2rad( [30.0 0.0 0.0] ) );
    T = [0.5 0.0 0.0];
    pcData1 = Loc2Glo( pcData0, R', T' );
    
    tmp = pcdownsample( pointCloud(pcData0'), 'gridAverage', CoarseRes );
    [SIdx, DD ]= knnsearch( pcData0', tmp.Location, 'k', 1 );
    RefData = pcData0(:, SIdx);
    MovData = pcData1(:, SIdx);
    
    params.DistThr = 0.002;
    params.MaxIter = Inf;
end
%%%%%%%%%%%% parallel computing to accelerate.
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)   % if exists no parallel computing handle, create one.
    poolsize = parpool;
else
    poolsize = p.NumWorkers;
end
%%%%%%%%%%%% start to recovery.
% tic;
DataLen = size(MovData, 2);
X = 1:1:DataLen;
C = [];
parfor i = 1 : 1 : (length(X)-1)
    k = [ repelem(i, length(X)-i); (i+1):1:length(X) ];
    C = [C k];
end
C = C';
% nCTime = toc;
% disp( sprintf('C generating Time = %.2fms', nCTime * 1000.0) );

% tic
EvalNum = [];
DistThr = params.DistThr;
Tf = zeros(3, 4, size(C, 1) );
MaxIter = min( [ params.MaxIter size(C, 1) ] ); 
parfor i = 1 : 1 : MaxIter
    [R, T, tmpDist] = RecoverRT( MovData, RefData, C(i, :) ); % RecoverRT( MovData( :, C(i, :) ), RefData( :, C(i, :) ) );
    EvalNum = [ EvalNum length( find( tmpDist < DistThr^2 ) )];
    Tf(:, :, i) = [R T];
end
% nPermuteTime = toc; 
% disp( sprintf('Permute Time = %.2fms', nPermuteTime * 1000.0) );
%%%%%%%%%%%%%%%% Refine the results.
% tic
[maxEvalNum Idx] = max(EvalNum);
R0 = Tf(:, 1:end-1, Idx);
T0 = Tf(:, end, Idx );
%%%%%%%%%%%%%%% Refine the transformations.
tmp = bsxfun( @plus, R0 * MovData, T0 ) - RefData;
tmpDist = tmp(1, :).^2 + tmp(2, :).^2 + tmp(3, :).^2;
Idx = find( tmpDist < DistThr^2 );
[R0, T0, tmpError ] = RecoverRT( MovData, RefData, Idx);
% nRefineTime = toc;
% disp( sprintf('Refine Time = %.2fms', nRefineTime * 1000.0) );
bTest = 1;
end

