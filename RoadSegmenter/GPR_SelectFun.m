function out = GPR_SelectFun(data, VIdx, tData, tDist, Hyp)
if nargin == 0
end
if size(data, 2) ~= 2
    data = data';
end
if size(data, 2) ~= 2
    error('GPRGrdSegFun_WithFixedHyp() is wrong!');
end
out = []; 
out = struct('data', [], 'ValidIdx', [], 'NIdx', [], 'hyp', [] );
out.data = data;
nLen = max(size(data)); 
NIdx = find( ~ismember([1:1:nLen], VIdx) ); 
for i = 1 : 1 : length(NIdx)
    id = NIdx(i);
    pt = data(id, :);
    xt = data(VIdx, 1); 
    yt = data(VIdx, 2); 
    [mu s2] = myGPRFun(Hyp, xt, yt, pt(1));
    Ratio = (mu - pt(2))./sqrt(s2) / tData;
    %%%%%%%%% update idx.
   if abs(Ratio) <= 1 & sqrt(s2) <= tDist 
       VIdx(end+1) = id; 
   end
end
NIdx = find( ~ismember([1:1:nLen], VIdx) ); 
out.ValidIdx = VIdx;
out.NIdx = NIdx;
out.hyp  = Hyp;
end