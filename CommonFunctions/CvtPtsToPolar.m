function [varargout]  = CvtPtsToPolar( pcData, varargin )
switch (nargin-1)
    case 0
        RadArray = [ 0.0 : 0.2 : 20.0 20.5 : 0.5 : 50.0 ];
        AngRes   = deg2rad(1.0);
    case 1
        AngRes = varargin{1};
        RadArray = [ 0.0 : 0.2 : 20.0 20.5 : 0.5 : 50.0 ];
    case 2
        RadArray = varargin{1};
        AngRes   = varargin{2};
    otherwise
        error('Invalid input!\n');
end
Angle = wrapTo2PiFun( atan2( pcData(2, :) , pcData(1, :) ) );
SegID = ceil( (Angle+1e-6)/ AngRes );
maxLen = ceil(2*pi/AngRes); 
SegID(SegID>maxLen) = maxLen;
USE_ANGLE_RANGE = 0; 
if USE_ANGLE_RANGE
    EffAngIdx = find( Angle <= deg2rad(90.0) | Angle >= deg2rad(270.0) ); 
end
Radius = sqrt( pcData(1, :).^2 + pcData(2, :).^2 );
BinID = -1 * ones(1, size(pcData, 2) );   % BinID = 1 stores invalid points, such as radius smaller than 4.5m or bigger than 5.0m.
for i = 1 : 1 : (length(RadArray)-1)
    Idx = find( Radius >= RadArray(i) & Radius < RadArray(i+1) );
    if USE_ANGLE_RANGE
        ind = ismember(Idx, EffAngIdx); 
        Idx = Idx(ind); 
    end
    BinID(Idx) = i;
end
if nargout == 2
    varargout{1} = SegID;
    varargout{2} = BinID;
end
if nargout == 3
    varargout{1} = SegID;
    varargout{2} = BinID;
    varargout{3} = find(BinID ~= -1);
end