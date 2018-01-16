function [varargout] = p2plICP(Md, MovData, RefData, Norm_Mov, Norm_Ref, Tf0, DistThr, AngThr )
if nargin == 0
    clc; close all;
    A = pcread('teapot.ply');  %pcread('E:\Dataset\Point Cloud\ply\Bunny\reconstruction\bun_zipper.ply');
    A = pcdownsample(A, 'gridAverage', 0.10 );
    MovData = double(A.Location');
    R0 = eul2rotm( deg2rad([20.0 10.0 10.0]) );
    T0 = [5 5 10.00 ]';
    RefData = Loc2Glo( MovData, R0', T0 );
    Sigma = 0.01;
    % rng(100);
    MovData = MovData + Sigma * randn(3,  size(MovData, 2) );
    NormalInfo_Ref = CalNormalsFun(RefData, 20, zeros(3, 1) );
    NormalInfo_Mov = CalNormalsFun(MovData, 20, zeros(3, 1) ); 
    Norm_Ref = cat(2, NormalInfo_Ref(:).normal);
    Norm_Mov = cat(2, NormalInfo_Mov(:).normal);
    Md = createns(RefData');
    Tf0 = [ eye(3) mean(RefData-MovData, 2); 0 0 0 1];
    DistThr = Inf;
    [NNIdx, DD] = knnsearch( RefData', RefData', 'k', 2);
    Res = mean(DD(:, 2));
    DistThr = Inf; 
    AngThr = cosd(45); 
    % ErrTol = [0.1 Res / 20];
    % RefData = RefData(:, 1:3000); 
end
bVerbose = false; % true; 
MaxIter = 20;
Tf = Tf0;
TotalTf.Tf = Tf;
ObjVal = [];
if isstruct(Norm_Mov)
     Norm_Mov = cat(2, Norm_Mov(:).normal);
end
if isstruct(Norm_Ref)
     Norm_Ref = cat(2, Norm_Ref(:).normal);
end
for Iter = 1 : 1 : MaxIter
    %%%%%%%%%%% establish correspondence.
    AftData = Loc2Glo( MovData, Tf(1:3, 1:3)', Tf(1:3, end) );
    [NNIdx, DD] = knnsearch( Md, AftData' );    
    Angle = sum( Norm_Ref(:, NNIdx) .* (Tf(1:3, 1:3) * Norm_Mov) );
    idx = find( abs(Angle) > AngThr & DD' < DistThr );
    s = AftData(:, idx );
    d = RefData(:, NNIdx(idx) );
    n = Norm_Ref( :, NNIdx(idx) );

    if Iter == 1
        ObjVal(end+1) = CalObj(Tf, s, d, n);
        if bVerbose 
            disp(sprintf('Iter = 00, ObjVal = %.3f', ObjVal) );
        end
    end
    %%%%%%%%%%% obtain rotation and translation via solving linear equations.
    tmpDiff = d - s;
    b = sum( n .* tmpDiff)';  
    a1 = n(3, :) .* s(2, :) - n(2, :) .* s(3, :);
    a2 = n(1, :) .* s(3, :) - n(3, :) .* s(1, :);
    a3 = n(2, :) .* s(1, :) - n(1, :) .* s(2, :);
    A = [a1' a2' a3' n'];
    %%%%%%%%%%% solve the problem || Ax - b ||, the solution also can be x = pinv( A' * A) * A' * b - x
    % x = pinv( A' * A) * A' * b;
    x = pinv(A) * b;
    
    [ dR, dT ] = ObtainTf(x);
    Tf(1:3, 1:3) = dR * Tf(1:3, 1:3);
    Tf(1:3, end) = dR * Tf(1:3, end) + dT;
    TotalTf(end+1).Tf = Tf;
    %%%%%%%%%%% check convergence condition.
    Err = max( norm(dR - eye(3)), norm(dT) );  % CalRT_Diff(TotalTf(end).Tf, TotalTf(end-1).Tf );
    if bVerbose
        ObjVal(end+1) = CalObj(Tf, s, d, n);
        tmpStr = sprintf('Iter = %02d, ObjVal = %.3f, Err =  %.6f', Iter, ObjVal(end), Err );
        disp(tmpStr);
    end
    if Err(1) <= 1e-6
        break;
    end
end
if nargout == 1
    varargout{1} = Tf;
end
if nargout == 2
    varargout{1} = Tf(1:3, 1:3); 
    varargout{2} = Tf(1:3, end); 
end
IS_SHOW = 1;
if IS_SHOW
    figure;
    hold on;
    view(3);
    showPointCloud(RefData', 'g');
    showPointCloud(MovData', 'r');
    AftData = Loc2Glo( MovData, Tf(1:3, 1:3)', Tf(1:3, end) );
    showPointCloud(AftData', 'b');
    
    if nargin == 0
        ErrR = RotationDiff( Tf(1:3, 1:3), R0 );
        ErrT = norm(Tf(1:3, end) - T0 );
        title(sprintf( 'Point-to-Plane ICP, ErrR = %.4fdegs ErrT = %.4f', rad2deg(ErrR), ErrT ));
    else
        title('Point-to-Plane ICP'); 
    end
end;
bTest = 1;
end

function varargout = ObtainTf(x)
Ang = x(3:-1:1);
R = eul2rotm(Ang');
T = x(4:6);
if nargout == 2
    varargout{1} = R;
    varargout{2} = T;
end
if nargout == 1
    varargout{1} = [R T; 0 0 0 1];
end
end

function [Err] = CalRT_Diff(Tf1, Tf0)
Err = zeros(2, 1);
Err(1) = RotationDiff(Tf1(1:3, 1:3), Tf0(1:3, 1:3));
Err(1) = rad2deg(Err(1));
Err(2) = norm( Tf1(1:3, end) - Tf0(1:3, end) );
end

function ObjVal = CalObj(M, s, d, n)
N = size(s, 2);
ObjVal = 0.0;
for i = 1 : 1 : N
    tmp = dot( M(1:3, 1:3)*s(:, i)+M(1:3,end) - d(:, i), n(:, i) );
    ObjVal = ObjVal + tmp^2;
end
end
