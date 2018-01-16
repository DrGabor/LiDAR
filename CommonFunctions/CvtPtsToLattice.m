function [ ImgInfo, varargout ] = CvtPtsToLattice( pcData, varargin )
switch (nargin-1)
    case 0
        Coverage = [-50.0 50.0; -50.0 50.0];
        Res = 0.2;
    case 1
        Coverage = [-50.0 50.0; -50.0 50.0];
        Res = varargin{1};
    case 2
        Coverage = varargin{1};
        Res = varargin{2};
    otherwise
        error('Invalid input!\n');
end
W_Wide = ceil( ( Coverage(2, 2) - Coverage(2, 1) ) / Res );
H_Wide = ceil( ( Coverage(1, 2) - Coverage(1, 1) ) / Res );    % width and height of the grid map image.
ImgInfo = [ H_Wide W_Wide];
if isempty(pcData)
    XI = [];
    YI = [];
else
    X = floor( pcData(1, :) / Res );
    Y = floor( pcData(2, :) / Res );
    maxX = Coverage(1, 2); 
    maxY = Coverage(2, 2); 
    XI = ceil(maxY/ Res) - Y;
    XI( XI < 1 | XI > W_Wide ) = -1;
    YI = ceil(maxX/ Res) - X;
    YI( YI < 1 | YI > H_Wide ) = -1;
end
if (nargout-1) == 1
    varargout{1} = [ XI; YI];
end
if (nargout-1) == 2
    varargout{1} = XI;
    varargout{2} = YI;
end
end

