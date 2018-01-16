function [varargout] = p2pSimpleICP(Md, MovData, RefData, Tf0, DistThr )
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
    Md = createns(RefData');
    Tf0 = [ eye(3) mean(RefData-MovData, 2); 0 0 0 1];
    DistThr = Inf;
    [NNIdx, DD] = knnsearch( RefData', RefData', 'k', 2);
    Res = mean(DD(:, 2));
    DistThr = Inf;
    AngThr = cosd(45);
end
bVerbose = false; % true;
nDim = size(RefData, 1); 
MaxIter = 30;
R = Tf0(1:nDim, 1:nDim);
T = Tf0(1:nDim, end); 
for Iter = 1 : 1 : MaxIter
    %%%%%%%%%%% establish correspondence.
    AftData = Loc2Glo( MovData, R', T );
    [NNIdx, DD] = knnsearch( Md, AftData' );
    idx = find( DD' < DistThr );
    MovIdx = idx;
    RefIdx = NNIdx(idx);
    %%%%%%%%%%% obtain rotation and translation via SVD.
    [ dR, dT ] = RegFun(RefData(:, RefIdx), AftData(:, MovIdx) );
    R = dR * R;
    T = dR * T + dT;
    %%%%%%%%%%% check convergence condition.
    Err = max( norm(dR - eye(nDim)), norm(dT) );  % CalRT_Diff(TotalTf(end).Tf, TotalTf(end-1).Tf );
%     str = sprintf( 'Iter = %02d, Err = %f\n', Iter, Err ); 
%     disp(str); 
    if Err(1) <= 1e-4
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
    grid on; 
    axis equal; 
    AftData = Loc2Glo( MovData, R', T );
    if nDim == 3
        view(3);
        showPointCloud(RefData', 'g');
        showPointCloud(MovData', 'r');
        showPointCloud(AftData', 'b');
    else
        plot(RefData(1, :), RefData(2, :), 'g.'); 
        plot(MovData(1, :), MovData(2, :), 'r.');
        plot(AftData(1, :), AftData(2, :), 'bo', 'markersize', 3);
    end
    
    if nargin == 0
        ErrR = RotationDiff( R, R0 );
        ErrT = norm(T - T0 );
        title(sprintf( 'Point-to-Point ICP, ErrR = %.4fdegs ErrT = %.4f', rad2deg(ErrR), ErrT ));
    else
        title('Point-to-Point ICP');
    end
end
bTest = 1;
end


