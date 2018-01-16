function [IDX, maxIdxArray] = Cluster(data, MinPts, epsilon, PtsNumThr, IS_SHOW)
if size(data, 2) ~= 2
    data = data';
end
if size(data, 2) ~= 2
    error('Input for Cluster() is wrong!');
end
IDX = [];
maxIdxArray = [];
if nargin == 4
    IS_SHOW = 0;
end
if size(data, 1) < PtsNumThr   %%%%%% too few points!
    return;
end
[ IDX, isNoise] =DBSCAN(data,epsilon,MinPts);
if IS_SHOW
    figure;
    hold on; 
    box on; 
    xlabel('X(meter)'); 
    ylabel('Y(meter)'); 
    PlotClusterinResult(data, IDX);
    title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);
end
ScanData = CalScanData(IDX);
if isempty(ScanData)
    IDX = []; 
    maxIdxArray = []; 
    return; 
end
[maxVal, idxOrder]= max(ScanData(2, :));
if maxVal < PtsNumThr
    return;
end
maxIdxArray = find(IDX == ScanData(1, idxOrder));
end

