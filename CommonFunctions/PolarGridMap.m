function PolarGM = PolarGridMap( pcData, varargin )
OpenParFor();
switch (nargin-1)
    case -1
        pcData = [ ( -20.0 : 0.1 : 40.0 ); 1.0 * ( -20.0 : 0.1 : 40.0 ) + 4.0; zeros( 1, length(( -20.0 : 0.1 : 40.0 )))];
        RadArray = [ 0.0 : 0.2 : 20.0 20.5 : 0.5 : 50.0 ];
        AngRes   = deg2rad(1.0);
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
if AngRes >= deg2rad(10.0)
    tmpStr = sprintf('The AngRes you define is %.2f degree, too big, is it the confusion of degree and radians? ', rad2deg(AngRes) );  
    warning(tmpStr); 
end
PolarGM = [];
[SegID BinID]  = CvtPtsToPolar( pcData, RadArray, AngRes );
%% construct polar grid map.
EffIdx = find( SegID > 0 & BinID > 0 );
SegID = SegID(EffIdx);
BinID = BinID(EffIdx);

SegNum = floor( 2 * pi / AngRes );
Idx = sub2ind( [SegNum length(RadArray)], SegID, BinID );
[C ia ic] = unique(Idx);
[r c] = ind2sub( [SegNum length(RadArray)], C );
GMIdx = [r; c];
TotalPixel = [];
s_elem = struct('Points', [], 'RawIdx', [], 'Gap', -Inf, 'MPts', [], 'ptLow', [] );
parfor i = 1 : 1 : length(C)
    x = GMIdx(1, i);
    y = GMIdx(2, :);
    tmp = s_elem;
    RealIdx = find( Idx == C(i) );
    RealIdx = EffIdx(RealIdx);
    tmp.RawIdx = RealIdx;
    Pts = pcData(:, RealIdx);
    tmp.Points = Pts;
    tmp.MPts = mean(Pts, 2);
    [~, minId] = min(Pts(3, :));
    tmp.ptLow = Pts(:, minId);
    if size(tmp.Points, 2) >= 2
        tmp.Gap = max(Pts(3, :)) - min(Pts(3, :));
    end
    TotalPixel = [TotalPixel tmp];
end
PolarGM = repmat( s_elem, SegNum, length(RadArray) );
for i = 1 : 1 : length(C)
    x = GMIdx(1, i);
    y = GMIdx(2, i);
    PolarGM(x, y) = TotalPixel(i);
end

%% Visualization part.
IS_SHOW = 0;
if IS_SHOW
    Gap = cat(1, PolarGM(:).Gap );
    Idx = find( Gap >= 0.5 );
    pcDataObs = cat(2, PolarGM(Idx).Points );
    figure;
    hold on; axis equal; grid on; view(3);
    showPointCloud(pointCloud(pcDataObs(1:3, :)') );
end
end

