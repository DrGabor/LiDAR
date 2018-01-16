function [Pts, varargout] = Loc2Loc(pcData, PoseRef, PoseAft )
if isequal(PoseRef, PoseAft)
    Pts = pcData;
    R0 = eye(3);
    R1 = R0;
    T0 = zeros(3, 1);
    T1 = T0;
else
    [Pts, R0, T0] = Loc2Glo( pcData, PoseRef);
    [Pts, R1, T1] = Glo2Loc( Pts, PoseAft);
end

R = R0 * R1;
T = inv(R1) * T0 + T1;
if (nargout-1) == 1
    varargout{1} = [R T];
end
if (nargout-1) == 2
    varargout{1} = R;
    varargout{2} = T;
end
end

