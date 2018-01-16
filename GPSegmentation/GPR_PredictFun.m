function [mu, s2, F, EffIdx, NffIdx] = GPR_PredictFun(hyp, data, xs, ys, tData, tDist)
if size(data, 2) ~= 2
    data = data';
end
if size(data, 2) ~= 2
    error('Input for GPR_PredictFun is wrong!');
end
if isrow(xs)
    xs = xs'; 
end
if isrow(ys)
    ys = ys'; 
end
x = data(:, 1);
y = data(:, 2);
[mu s2] = myGPRFun(hyp, x, y, xs);
F = [mu+tData*sqrt(s2); flipdim(mu-tData*sqrt(s2),1)];
if ~isempty(ys)
    Ratio = (mu - ys)./sqrt(s2) / tData;
    %%%%%%%%% update idx.
    EffIdx = find( abs(Ratio) <= 1 & sqrt(s2) <= tDist );
    NffIdx = find( ~ismember( 1:1:length(xs), EffIdx) );
    if isrow(NffIdx)
        NffIdx = NffIdx';
    end
else
    EffIdx = [];
    NffIdx = [];
end
end

