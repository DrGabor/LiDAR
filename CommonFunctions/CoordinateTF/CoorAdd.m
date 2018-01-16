%% Coordinate add, the pose0 is global coordinate transform, after addiing IncrePose, obtaining the global coordinate transform Pose1.
function [varargout] = CoorAdd( Pose0, IncrePose )
[R0, T0] = ExtractRT( Pose0 );
[dR, dT] = ExtractRT( IncrePose );
% supose we have points in Pose0 and Pose1 as Pts0 and Pts1, Pts0 adds dR
% and dT obtaining Pts1. Pts1 = inv(dR) * Pts0 + dT = R1 * ( inv(R0) * Pts0 + T0 - T1 ), we can obtain:   
R1 = dR' * R0;
T1 = T0 - R1' * dT;
% T1 = R0'* dT + T0;
if nargout == 1
    varargout{1} = [R1 T1];
end
if nargout == 2 
    varargout{1} = R1;
    varargout{2} = T1;
end
end

