function [GrdIdx, ObsIdx, UnkownIdx, GrdPts] = GPSegFun(pcData, params)
if nargin == 0
    clc; close all;
    DataFolder = 'D:\Data\Sudoku\Record-2017-10-26-09-03-48(Match2016)\BinaryData';
    nFrm = 655;
    DataDir = fullfile(DataFolder, sprintf('Binary%d.txt', nFrm));
    pcData = HDLAnalyserNew(DataDir);
    pcData = pcData(1:3, :);
    params = []; 
    params.GapThr = [0.10 1.0]; 
    params.RadArray = 0.0 : 2.0 : 80.0; 
    params.AngRes   = deg2rad(5.0); 
    params.GroundH  = -1.8; 
    params.IS_SHOW = 1; 
end
if nargin == 1
    params = []; 
    params.GapThr = [0.10 1.0]; 
    params.RadArray = 0.0 : 2.0 : 80.0; 
    params.AngRes   = deg2rad(5.0); 
    params.GroundH  = -1.8; 
    params.IS_SHOW = 0; 
end
[SegID, BinID, EffIdx] = CvtPtsToPolar(pcData, params.RadArray, params.AngRes);
SegID = SegID(EffIdx);
BinID = BinID(EffIdx);
data = pcData(:, EffIdx);
maxAngLen = ceil(2*pi/AngRes);
maxBinLen = length(RadArray);
ind = sub2ind( [maxAngLen maxBinLen], SegID, BinID );
[B I] = sort(ind);
[C ia ic] = unique(B);
ss = struct('Points', [], 'RawIdx', [], 'Gap', -Inf, 'MPts', -Inf(3, 1), 'CoorI', [] );
TotalPixel = [];
[seg bin] = ind2sub([maxAngLen maxBinLen], C);
parfor i = 1 : 1 : length(C)   % C = ind(ia) or
    id0 = ia(i);
    if i == length(C)
        id1 = length(ind);
    else
        id1 = ia(i+1) - 1;
    end
    Idx = I(id0:1:id1);
    if length(Idx) > 10
        bTest = 1;
    end
    a = unique(ind(Idx));
    if length(a) ~= 1
        error('wrong!');
    end
    tmp = ss;
    tmp.RawIdx = EffIdx(Idx);
    tmp.CoorI = [seg(i); bin(i)];
    TotalPixel = [TotalPixel tmp];
end
CoorI = cat(2, TotalPixel(:).CoorI);
ObsIdx = [];
UnkownSeg0 = [];
ss = struct('data', [] );
RMInfo = repmat(ss, maxAngLen, maxBinLen);
verbose = 0;
UnkownIdx = [];
rand('state',0);
for AngID = 1 : 1 : maxAngLen
    SelIdx = find(CoorI(1, :) == AngID );
    SelPts = [];
    for i = 1 : 1 : length(SelIdx)
        tmp = TotalPixel(SelIdx(i));
        Pts = pcData(:, tmp.RawIdx);
        [~, id] = min(Pts(3, :));
        SelPts(:, end+1) = Pts(:, id);
        bTest = 1;
    end
    EIdx = cat(2, TotalPixel(SelIdx).RawIdx);
    if size(SelPts, 2) < 10
        UnkownSeg0 = [UnkownSeg0 AngID];
        str = sprintf('AngID = %03d/%03d, too few points, continue!', AngID, maxAngLen);
        if verbose
            disp(str);
        end
        continue;
    end
    %%%%%%%
    Dist = sqrt( SelPts(1, :).^2 + SelPts(2, :).^2 );
    Pts = [Dist; SelPts(3, :)];
    DistThr = 0.2;
    [V, L, inliers] = ransacfit2Dline(Pts, 0.1, 0);
    Dist = ( Pts(1, :)*V(1) + Pts(2, :)*V(2)+V(3) ) / norm(V(1:2));
    %%%%% y = kx + b.
    k = -V(1)/V(2);
    b = -V(3)/V(2);
    testRad = 5.0;
    H0 = k*testRad + b;
    idx = find(abs(Dist) <= DistThr);
    if length(idx) < 10 | ( H0 > (params.GroundH + 0.3) | H0 < (params.GroundH - 0.3) )  % the orignal ground height is -1.8m.
        UnkownSeg0 = [UnkownSeg0 AngID];
        if length(idx) < 10
            str = sprintf('AngID = %03d/%03d, too few initial points, continue!', AngID, maxAngLen);
        else
            str = sprintf('AngID = %03d/%03d, ransac finds a wrong line, continue!', AngID, maxAngLen);
        end
        if verbose
            disp(str);
        end
        continue;
    end
    str = sprintf('AngID = %03d/%03d', AngID, maxAngLen);
    if verbose
        disp(str);
    end
    tData = 2.0;
    tDist = 0.15;
    % save('./GInfo.mat');
    % out = GPRGrdSegFun(Pts, idx', tData, tDist, [], 0);
    Hyp.cov = log([6.0 sqrt(1.3298)]);
    Hyp.mean = [0 0];
    Hyp.lik = log(0.1);
    out = GPR_SelectFun(Pts, idx, tData, tDist, Hyp);
    VIdx = out.ValidIdx;
    xs = Pts(1, :);
    ys = Pts(2, :);
    data_t = out.data(out.ValidIdx, :);
    [mu, s2, F, EffIdx, NffIdx] = GPR_PredictFun(out.hyp, data_t, xs, ys, tData, tDist);
    [H0, ~] = myGPRFun(out.hyp, data_t(:, 1), data_t(:, 2), testRad);
    if ( H0 > -1.5 | H0 < -2.1 )
        UnkownSeg0 = [UnkownSeg0 AngID];
        continue;
    end
    CorrectH = ys';
    CorrectH(NffIdx) = mu(NffIdx);
    for i = 1 : 1 : length(SelIdx)
        tmp = TotalPixel(SelIdx(i));
        Pts = pcData(:, tmp.RawIdx);
        tmpDiff = Pts(3, :) - CorrectH(i);
        idx = find(tmpDiff >= GapThr(1) & tmpDiff <= GapThr(2) );
        ObsIdx = [ObsIdx tmp.RawIdx(idx)];
        tmpIdx = find(tmpDiff > GapThr(2));
        UnkownIdx = [UnkownIdx tmp.RawIdx(tmpIdx)];
        bTest = 1;
        if size(Pts, 2) > 1
            MPt = mean(Pts')';
        else
            MPt = Pts;
        end
        MPt(3) = CorrectH(i);
        RMInfo( tmp.CoorI(1), tmp.CoorI(2) ).data = MPt;
    end
    bTest = 1;
end
%% check the unkown segments.
GrdPts = cat(2, RMInfo(:).data);
GrdMd = createns(GrdPts');
UnkownSeg = [];
IS_SHOW_DETAIL = 0;
NewObsIdx = [];
rArray = [2.0 5.0 10.0 15.0 20.0 25.0 30.0];
for id = 1:1:length(UnkownSeg0)
    AngID = UnkownSeg0(id);
    SelIdx = find(CoorI(1, :) == AngID );
    TTIdx = cat(2, TotalPixel(SelIdx).RawIdx);
    SelPts = [];
    for i = 1 : 1 : length(SelIdx)
        tmp = TotalPixel(SelIdx(i));
        Pts = pcData(:, tmp.RawIdx);
        if size(Pts, 2) > 1
            MPt = mean(Pts')';
        else
            MPt = Pts;
        end
        for rId = 1:1:length(rArray)
            [NNIdx, DD] = rangesearch(GrdMd, MPt', rArray(rId));
            if ~isempty(NNIdx{1})
                break;
            end
        end
        if isempty(NNIdx{1})
            UnkownIdx = [UnkownIdx tmp.RawIdx];
            UnkownSeg = [UnkownSeg AngID];
            continue;
        end
        NNIdx = NNIdx{1};
        tmpPts = GrdPts(:, NNIdx);
        H = mean(tmpPts(3, :));
        tmpDiff = Pts(3, :) - H;
        % tmpDiff = abs(tmpDiff);
        idx = find(tmpDiff >= GapThr(1) & tmpDiff <= GapThr(2) );
        NewObsIdx = [NewObsIdx tmp.RawIdx(idx)];
        tmpIdx = find(tmpDiff > GapThr(2));
        UnkownIdx = [UnkownIdx tmp.RawIdx(tmpIdx)];
    end
    if IS_SHOW_DETAIL
        figure;
        hold on;
        grid on;
        axis equal;
        view(3);
        data = pcData(:, NewObsIdx);
        showPointCloud(data', 'r');
        OIdx = find(~ismember(TTIdx, NewObsIdx));
        data = pcData(:, TTIdx(OIdx));
        showPointCloud(data', 'g');
    end
    bTest = 1;
end
ObsIdx = [ObsIdx NewObsIdx];
Dist = sqrt(pcData(1, :).^2 + pcData(2, :).^2);
NIdx = find(Dist > RadArray(end) );
GrdIdx = find(~ismember(1:1:size(pcData,2), ObsIdx) & ~ismember(1:1:size(pcData,2), UnkownIdx) & ~ismember(1:1:size(pcData,2), NIdx) );
GrdMd = createns(GrdPts(1:2, :)');
GrdAddPts = []; 
for i = 1 : 1 : size(RMInfo, 1)
    x = i*ones(1, size(RMInfo, 2));
    y = 1:1:size(RMInfo, 2);
    ind = sub2ind(size(RMInfo), x, y);
    Tmp = cat(2, RMInfo(ind).data);
    if isempty(Tmp)
        SelIdx = find(CoorI(1, :) == i );
        Tmp = [];
        for i = 1 : 1 : length(SelIdx)
            tmpInfo = TotalPixel(SelIdx(i));
            coorI = tmpInfo.CoorI;
            angId = coorI(1);
            radId = coorI(2);
            ang = (angId - 0.5)*AngRes;
            r = RadArray(radId) + 0.5*RadRes;
            pt = [r*cos(ang); r*sin(ang)];
            
            for rId = 1:1:length(rArray)
                [NNIdx, DD] = rangesearch(GrdMd, pt', rArray(rId));
                if ~isempty(NNIdx{1})
                    break;
                end
            end
            NNIdx = NNIdx{1};
            tmpPts = GrdPts(:, NNIdx);
            h = mean(tmpPts(3, :));
            pt(end+1) = h;
            %                 Pts = pcData(:, tmpInfo.RawIdx);
            %                 [~, id] = min(Pts(3, :));
            GrdAddPts(:, end+1) = pt;
            bTest = 1;
        end
    end
end
GrdPts = [GrdPts GrdAddPts];
if params.IS_SHOW
    figure;
    hold on;
    grid on;
    view(3);
    
    if ~isempty(GrdIdx)
        showPointCloud(pcData(1:3, GrdIdx)', 'g');
    end
    if ~isempty(UnkownIdx)
        showPointCloud(pcData(1:3, UnkownIdx)', 'k');
    end
    if ~isempty(ObsIdx)
        showPointCloud(pcData(1:3, ObsIdx)', 'r');
    end
    title('Segmentation via GPR');
    xlabel('X/m');
    ylabel('Y/m');
    zlabel('Z/m');
    figure; 
    hold on; 
    grid on; 
    axis equal; 
    pcshow(pcData(1:3, :)'); 
    plot3(GrdPts(1, :), GrdPts(2, :), GrdPts(3, :), 'ro', 'markersize', 3, 'markerfacecolor', 'r', ...
        'markeredgecolor', 'y'); 
    xlabel('X/m');
    ylabel('Y/m');
    zlabel('Z/m');
    title('Road Model'); 
end
end
