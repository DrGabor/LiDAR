%% pair-wise match to build a cohenent map.
function RoPSInfo = CalRoPSFun(pcDataRaw, ResArray)
%%%%%%%%%%% Calculate RoPS Features.
if isempty(ResArray)
    Res0 = CalPointsRes(pcDataRaw);
    FineRes   =  1.0 * Res0;
    CoarseRes = 10.0 * Res0;
else
    FineRes   = ResArray(1);
    CoarseRes = ResArray(2);
end
%%%%%%%%%%%%%%%%% Fine sampling the raw points.
FineCloud0 = pcdownsample( pointCloud(pcDataRaw'), 'gridAverage', FineRes );
pcData0 = FineCloud0.Location';

%%%%%%%%%%%%%%%%% Coarse sampling to generate candidate points and also
%%%%%%%%%%%%%%%%% calculate reliable points.
CoarseCloud0 = pcdownsample( pointCloud(pcDataRaw'), 'gridAverage', CoarseRes );
[SIdx0, DD ]= knnsearch( pcData0', CoarseCloud0.Location, 'k', 1 );
Tao = 1.05;
[~, ~, ~, EigArrays] = pcNormals(pcData0 );
tmpIdx = find( EigArrays(3, SIdx0) ./ EigArrays(2, SIdx0) >= Tao );
SIdx0 = SIdx0(tmpIdx);

%%%%%%%%%%%%%%%%% use RPoS features.
pointcloud0 = RPoS(pcData0, SIdx0 );

%%%%%%%%%%%%%%%% assign into structure.
RoPSInfo.FineRes = FineRes;
RoPSInfo.CoarseRes = CoarseRes;
RoPSInfo.pcDataRaw = pcDataRaw;
RoPSInfo.vertices = pcData0;
RoPSInfo.keyIdx = SIdx0;
RoPSInfo.RoPS = cat( 2, pointcloud0.RoPS{:} )';
RoPSInfo.LRF = pointcloud0.LRF;
end
