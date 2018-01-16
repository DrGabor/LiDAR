function [Pts, varargout] = Glo2Loc( pcData, varargin )
[ Dim, DataLen] = size(pcData);
if Dim ~= 2 && Dim ~= 3
    error('the input pcData must be 2 * N or 3 * N format!');
end
if (nargin - 1) == 2
    R = varargin{1};
    T = varargin{2};
end
if (nargin - 1) == 1
    Pose = varargin{1};
    [R T] = ExtractRT(Pose);
end
Pts = R * ( bsxfun( @minus, pcData, T ) );
if (nargout-1) == 1
    varargout{1} = [inv(R) -R*T];
end
if (nargout-1) == 2
    varargout{1} = inv(R);
    varargout{2} = -R*T;
end
end

