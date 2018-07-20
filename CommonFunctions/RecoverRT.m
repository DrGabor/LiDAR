%% data must be the Dim * N format.
function [R ,T, Errors] = RecoverRT(MovPts, RefPts, varargin)
if nargin == 0
    clc;
    MovPts = rand(3, 5);
    R0 = eul2rotm( deg2rad([30.0 2.5 45.0]) );
    T0 = [ 0.5; 0.2; 0.0];
    RefPts = bsxfun(@plus, R0*MovPts, T0);
    SelIdx = 1 : 1 : size(RefPts, 2);
end
%%%%%%%%%%%% check inputs' validation.
if size(MovPts, 1) ~= 2 && size(MovPts, 1) ~= 3
    error('MovPts must be the Dim * N format.');
end
if size(RefPts, 1) ~= 2 && size(RefPts, 1) ~= 3
    error('RefPts must be the Dim * N format.');
end
if size(RefPts, 2) ~= size(MovPts, 2)
    error('RefPts and MovPts should be same size.');
end

if (nargin-2) == 0
    SelIdx = 1 : 1 : size(RefPts, 2);
end
if (nargin-2) == 1
    SelIdx = varargin{1};
end
Dim = size(MovPts, 1);
MovPts0 = MovPts;
RefPts0 = RefPts;
MovPts = MovPts(:, SelIdx);
RefPts = RefPts(:, SelIdx);
%% Using SVD to recovery the transformation: RefPts = R * MovPts + T.
mmRef = mean(RefPts, 2);
mmMov = mean(MovPts, 2);
MovPts = bsxfun(@minus, MovPts, mmMov );
RefPts = bsxfun(@minus, RefPts, mmRef );
K = MovPts*RefPts';
[U, S, V] = svd(K);
M = eye(Dim);
M(end,end) = det(V*U');
R = V*M*U';
T = mmRef - R*mmMov;
tmp = bsxfun( @plus, R * MovPts0, T ) - RefPts0;
if Dim == 2
    Errors = tmp(1, :).^2 + tmp(2, :).^2;
else
    Errors = tmp(1, :).^2 + tmp(2, :).^2 + tmp(3, :).^2;
end
bTest = 1;
end





