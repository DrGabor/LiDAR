% the input pcData must be Dim * N format.
function [Pts, varargout] = Loc2Glo( pcData, varargin )
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
    Pts = bsxfun( @plus, inv(R) * pcData, T );
    if (nargout-1) == 1
        varargout{1} = [R T];
    end
    if (nargout-1) == 2
        varargout{1} = R;
        varargout{2} = T;
    end
end

