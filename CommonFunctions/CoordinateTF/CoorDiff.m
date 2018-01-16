%% Coordinate Difference, the pose0 and pose1 are global coordinate transform, calculate a coordinate transform which transforms pose0 to pose1.
function [varargout] = CoorDiff(Pose1, Pose0)
[R0, T0] = ExtractRT( Pose0 );
[R1, T1] = ExtractRT( Pose1 );
% suppose the points in Pose1 is Pts1, points in Pose0 is Pts0; the dR and dT is what we wanted: the points in Pose0 adds dR and dT, obtaining the points in Pose1.
% Use above relationship, we can get, Pts1 = inv(dR) * Pts0 + dT = R1( inv(R0) * Pts0 + T0 - T1 ), we can get following:  
dR = R0 * R1';
dT = R1 * (T0-T1);
if nargout == 1 
    varargout{1} = [dR dT];
end
if nargout == 2 
    varargout{1} = dR;
    varargout{2} = dT;
end
end

