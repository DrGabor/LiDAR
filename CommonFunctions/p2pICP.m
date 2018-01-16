function [varargout] = p2pICP(Md, MovData, RefData, Norm_Mov, Norm_Ref, Tf0, DistThr, AngThr )
if nargin == 0
    clc; close all;
    A = pcread('teapot.ply');  %pcread('E:\Dataset\Point Cloud\ply\Bunny\reconstruction\bun_zipper.ply');
    A = pcdownsample(A, 'gridAverage', 0.10 );
    MovData = double(A.Location');
    R0 = eul2rotm( deg2rad([20.0 10.0 10.0]) );
    T0 = [5 5 10.00 ]';
    RefData = Loc2Glo( MovData, R0', T0 );
    Sigma = 0.01;
    rng(100);
    MovData = MovData + Sigma * randn(3,  size(MovData, 2) );
    [Norm_Ref] = CalNormalsFun(RefData, 20, zeros(3, 1) );
    [Norm_Mov] = CalNormalsFun(MovData, 20, zeros(3, 1) );
    Md = createns(RefData');
    Tf0 = [ eye(3) mean(RefData-MovData, 2); 0 0 0 1];
    DistThr = Inf;
    [NNIdx, DD] = knnsearch( RefData', RefData', 'k', 2);
    Res = mean(DD(:, 2));
    DistThr = Inf;
    AngThr = cosd(45);
end
bVerbose = false; % true;
MaxIter = 20;
Dim = size(RefData, 1); 
R = Tf0(1:Dim, 1:Dim);
T = Tf0(1:Dim, end); 
if isstruct(Norm_Mov)
     Norm_Mov = cat(2, Norm_Mov(:).normal);
end
if isstruct(Norm_Ref)
     Norm_Ref = cat(2, Norm_Ref(:).normal);
end
for Iter = 1 : 1 : MaxIter
    % Iter
    %%%%%%%%%%% establish correspondence.
    AftData = Loc2Glo( MovData, R', T );
    [NNIdx, DD] = knnsearch( Md, AftData' );
    if ~isempty(Norm_Ref) && ~isempty(Norm_Mov) && ~isempty(AngThr)
        Angle = sum( Norm_Ref(:, NNIdx) .* (R * Norm_Mov) );
        idx = find( abs(Angle) > AngThr & DD' < DistThr );
    else
        idx = find(DD' < DistThr); 
    end
    MovIdx = idx;
    RefIdx = NNIdx(idx);
    %%%%%%%%%%% obtain rotation and translation via SVD.
    [ dR, dT ] = RegFun(RefData(:, RefIdx), AftData(:, MovIdx) );
    R = dR * R;
    T = dR * T + dT;
    %%%%%%%%%%% check convergence condition.
    Err = max( norm(dR - eye(Dim)), norm(dT) );  % CalRT_Diff(TotalTf(end).Tf, TotalTf(end-1).Tf );
%     str = sprintf( 'Iter = %02d, Err = %f\n', Iter, Err ); 
%     disp(str); 
    if Err(1) <= 1e-6
        break;
    end
end
if nargout == 1
    varargout{1} = [R T];
end
if nargout == 2
    varargout{1} = R;
    varargout{2} = T;
end
IS_SHOW = 0;
if IS_SHOW
    figure;
    hold on;
    view(3);
    showPointCloud(RefData', 'g');
    showPointCloud(MovData', 'r');
    AftData = Loc2Glo( MovData, R', T );
    showPointCloud(AftData', 'b');
    
    if nargin == 0
        ErrR = RotationDiff( R, R0 );
        ErrT = norm(T - T0 );
        title(sprintf( 'Point-to-Point ICP, ErrR = %.4fdegs ErrT = %.4f', rad2deg(ErrR), ErrT ));
    else
        title('Point-to-Point ICP');
    end
end;
bTest = 1;
end


