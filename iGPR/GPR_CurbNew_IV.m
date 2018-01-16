function [out] = GPR_CurbNew_IV(data, ValidIdx, options0, InitHyp)
PtsNumThr = options0.PtsNumThr;
epsilon = options0.epsilon;
MinPts = options0.MinPts;
tData = options0.tData; 
tDist = options0.tDist; 
IS_SHOW = options0.IS_SHOW; 
if size(data, 2) ~= 2
    data = data';
end
if size(data, 2) ~= 2
    error('Input for GPR_CurbNew_IV is wrong!');
end
out = struct('data', [], 'ValidIdx', [], 'NIdx', [], 'hyp', [], 'nIter', [], 'basisVector', [] );
out.data = data;
if isempty(data)
    error('data is empty!');
end
basisVector = [];
%%%%%%%%%% DBSCAN clustering
XOffset = 0.0;
data = bsxfun(@plus, data, [XOffset 0.0]);
if isempty(ValidIdx)
    [~, maxIdxArray] = Cluster(data, MinPts, epsilon, PtsNumThr);
    if isempty(maxIdxArray)
        error('cluster failed! point number is smaller than PtsNumThr!');
    end
    ValidIdx = maxIdxArray;
    % IS_SHOW = 0;
end
NCheckIdx = find(~ismember(1:1:size(data, 1), ValidIdx) )';
%%%%%%%%%%%%%% start Gaussian Process Regression.
if isempty(InitHyp)
    InitHyp.mean = [0; 0];
    InitHyp.cov  = [0; 0];
    InitHyp.lik  = log(0.1);
end
% covfunc = {'covSum', {@covSEiso, @covNNone} }; covSEiso
covfunc = {@covSEiso};
hyp_cov = InitHyp.cov;
likfunc = @likGauss;  % Gaussian likelihood
hyp_lik = InitHyp.lik;
USE_ZERO_MEAN = 0;
if USE_ZERO_MEAN
    meanfunc = [];
    hyp_mean = [];
else
    nPoly = 1;
    if nPoly == 0
        meanfunc = @meanConst; % @meanConst;
        hyp_mean = 0.0;
    else
        mean_poly = {@meanPoly, nPoly};
        hypp = zeros(nPoly, 1);
        meanfunc = {'meanSum', {mean_poly, @meanConst}}; % @meanConst;
        hyp_mean = InitHyp.mean; % [hypp; 0.0];
    end
end
HypArray = struct('mean', hyp_mean, 'cov', hyp_cov, 'lik', hyp_lik);
nCount = 0;
SampleRate = 5;
USE_MY_GRADIENT = 1;
if USE_MY_GRADIENT
    options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region', ... % 'trust-region'
        'display', 'off', 'MaxIter', 100, 'TolFun', 1e-3, 'TolX', 1e-2 );
end
USE_MY_GPR = 1;
for id = 1 : 1 : 40
    x = data(ValidIdx, 1);
    y = data(ValidIdx, 2);
    if length(x) < 100
        CenPts = [x y]';
    else
        CenPts = ClusterRawFun(x, y);     %%%%%% trainning data, i.e. support vector.
    end
    hyp0 = HypArray(end);
    if USE_MY_GRADIENT
        X0 = [hyp0.mean; hyp0.cov];
        sn = 0.1;
        if id == 4
            bTest = 1;
        end
        if size(CenPts, 2) < 5
            continue;
        end
        xEst = fminunc(@myFunNew, X0, options, CenPts(1, :)', CenPts(2, :)', sn);
        hyp = hyp0;
        hyp.mean = xEst(1:2);
        hyp.cov  = xEst(3:4);
        tmp = [];
        tmp.vec = CenPts;
        basisVector = [basisVector tmp];
    else
        hyp = minimize(hyp0, @gp, max(1, ceil(30*0.9^(id-1))), @infGaussLik, meanfunc, covfunc, likfunc, CenPts(1, :)', CenPts(2, :)');
        hyp.lik = max(hyp_lik, hyp.lik);
    end
    HypArray = [HypArray hyp];
    trainData = [x(1:SampleRate:end) y(1:SampleRate:end)];
    %%%%%%% delete out noise data in valid index, this is computational intensive.
    if USE_MY_GPR
        [mu, s2, F, EffIdx, NffIdx] = GPR_PredictFun(hyp, trainData, x, y, tData, tDist);
    else
        [mu s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x(1:SampleRate:end), y(1:SampleRate:end), x);
    end
    %%%%%%%%% update idx.
    ValidIdx = ValidIdx(EffIdx);   % update valid idx.
    NCheckIdx = find(~ismember(1:1:size(data, 1), ValidIdx))';   % these index need to be further checked.
    
    %%%%%%% find out new valid data in NCheckIdx.
    if isempty(ValidIdx)
        continue;
    end
    trainData = [data(ValidIdx(1:SampleRate:end), 1) data(ValidIdx(1:SampleRate:end), 2)];
    if isempty(NCheckIdx)
        continue;
    end
    xs = data(NCheckIdx, 1);
    ys = data(NCheckIdx, 2);
    if USE_MY_GPR
        [mu, s2, F, EffIdx, ~] = GPR_PredictFun(hyp, trainData, xs, ys, tData, tDist);
    else
        [mu s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y, xs );
    end
    NewValidIdx = NCheckIdx(EffIdx);
    
    %%%%%%%%%%% calculate envelop and plot
    if IS_SHOW
        %%%%%%%%%%%%% calculate envelop.
        Xs = ( (min(x)-10.0):0.1:(max(x) + 10.0) )';
        if USE_MY_GPR
            [Mu, S2, F, ~, ~] = GPR_PredictFun(hyp, trainData, Xs, [], tData, tDist);
        else
            [Mu S2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y, Xs);
        end
        figure;
        hold on;
        grid on;
        axis equal;
        %%%%%%%%%% draw envelop.
        fill([Xs; flipdim(Xs,1)], F, [7 7 7]/8)
        plot(Xs, Mu);
        xlabel('X(meter)'); 
        ylabel('Y(meter)'); 
        box on; 
    end
    %%%%%%%%%% update valid idx.
    OldValidIdx = ValidIdx;
    ValidIdx = [ValidIdx; NewValidIdx];
    NCheckIdx = find(~ismember(1:1:size(data, 1), ValidIdx))';
    str = sprintf('Iter = %02d, ValidIdx = %03d, NCheckIdx = %03d, AddedIdx = %03d, DeleteIdx = %03d', ...
        id, length(ValidIdx), length(NCheckIdx), length(NewValidIdx), length(NffIdx));
    if IS_SHOW
        plot(data(OldValidIdx, 1), data(OldValidIdx, 2), 'gh', 'MarkerSize', 10 );
        plot(data(NewValidIdx, 1), data(NewValidIdx, 2), 'ks', 'MarkerSize', 10 );
        plot(data(NCheckIdx, 1), data(NCheckIdx, 2), 'rx', 'MarkerSize', 10 );
        title(str);
        disp(str);
    end
    
    if length(NewValidIdx) <= 0 & length(NffIdx) <= 1
        if IS_SHOW
            disp('Converge!');
        end
        nCount = nCount + 1;
    end
    if nCount >= 3
        break;
    end
    bTest = 1;
end
Xs = ( (min(x)-0.0):0.1:(max(x) + 10.0) )';
if USE_MY_GPR
    [Mu S2] = myGPRFun(hyp, x, y, Xs);
else
    [Mu S2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y, Xs);
end
F = [Mu+tData*sqrt(S2); flipdim(Mu-tData*sqrt(S2),1)];
Xs = Xs - XOffset;
meanData = cat(2, HypArray(2:end).mean);
covData  = exp(cat(2, HypArray(2:end).cov));
likData  = exp(cat(2, HypArray(2:end).lik) );
out.ValidIdx = ValidIdx;
out.NIdx = NCheckIdx;
out.hyp  = HypArray(end);
out.nIter = length(HypArray);
out.basisVector = basisVector;
if IS_SHOW
    figure;
    hold on;
    grid on;
    axis equal;
    box on; 
    if ~isempty(meanData)
        p0 = meanData(1:end-1, end);
        p1 = meanData(end, end);
        p = [p0(end:-1:1); p1];
        y = polyval(p, Xs);
        plot(Xs, y, 'k' );
    end
    plot(data(ValidIdx, 1), data(ValidIdx, 2), 'b.');
    plot(data(NCheckIdx, 1), data(NCheckIdx, 2), 'r.');
    plot(Xs, Mu, 'g');
    title('Filtered data');
end

%%%%%%%%%% draw hyperparameters.
if IS_SHOW
    if ~isempty(meanData)
        figure;
        hold on;
        for i = 1 : 1 : size(meanData, 1)
            data = meanData(i, :);
            subplot(size(meanData, 1), 1, i);
            hold on;
            grid on;
            plot(data, 'b*--');
            str = sprintf( 'Hyperparameter of mean, id = %02d', i);
            title(str);
        end
    end
    if ~isempty(covData)
        figure;
        hold on;
        for i = 1 : 1 : size(covData, 1)
            data = covData(i, :);
            subplot(size(covData, 1), 1, i);
            hold on;
            grid on;
            plot(data, 'b*--');
            str = sprintf( 'Hyperparameter of cov, id = %02d', i);
            title(str);
        end
    end
    if ~isempty(likData)
        figure;
        hold on;
        grid on;
        plot(likData, 'b*--');
        title('Hyperparameter of lik');
    end
end
bTest = 1;
end
function CenPts = ClusterRawFun(x, y)
if ~iscolumn(x)
    x = x';
end
if ~iscolumn(y)
    y = y';
end
IS_SHOW = 0;
if nargin == 0
    IS_SHOW = 1;
end
[~, idx] = sort(x);
data = [x(idx)'; y(idx)'];
tmpDiff = data(:, 2:end) - data(:, 1:end-1);
Dist = sqrt(tmpDiff(1, :).^2 + tmpDiff(2, :).^2);
CumDist = [1e-3 cumsum(Dist)];
if CumDist(end) > 50
    DistGap = 2.0;
else
    DistGap = 1.0;
end
tmp = ceil(CumDist/DistGap);
[~, b] = unique(tmp);
CenPts = [];
for i = 1 : 1 : (length(b)-1)
    if i == 57
        bTest = 1;
    end
    idx0 = b(i);
    if i == length(b)
        idx1 = size(data, 2);
    else
        idx1 = b(i+1) - 1;
    end
    Pts = data(:, idx0:1:idx1)';
    if size(Pts, 1) == 1
        pt = Pts;
    else
        pt = mean(Pts);
    end
    CenPts(:, end+1) = pt';
    bTest = 1;
end

if IS_SHOW
    figure;
    hold on;
    grid on;
    axis equal;
    plot(data(1, :), data(2, :), 'b.--');
    plot(data(1, 1), data(2, 1), 'rp' );
    plot(data(1, end), data(2, end), 'kp' );
    plot(CenPts(1, :), CenPts(2, :), 'mo', 'MarkerSize', 10);
    for i = 1 : 1 : size(CenPts, 2)
        text(CenPts(1, i), CenPts(2, i), num2str(i, '%02d'));
    end
    str = sprintf('DistGap = %.1f, CenPts = %02d', DistGap, size(CenPts, 2) );
    title(str);
end
end