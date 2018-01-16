function [varargout] = CvtLatticeToPts(varargin)
if nargin == 3
    XI = varargin{1}(1, :);
    YI = varargin{1}(2, :);
    Coverage = varargin{2};
    Res = varargin{3};
end
if nargin == 4
    XI = varargin{1};
    YI = varargin{2};
    Coverage = varargin{3};
    Res = varargin{4};
end
maxX = Coverage(1, 2);
maxY = Coverage(2, 2);
y = Res*(ceil(maxY/Res) - XI);
x = Res*(ceil(maxX/Res) - YI);
if nargout == 2
    varargout{1} = x;
    varargout{2} = y;
end
if nargout == 1
    varargout{1} = [x; y];
end

end