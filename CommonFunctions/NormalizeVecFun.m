function [Vec1] = NormalizeVecFun(Vec0)
n = size(Vec0, 2); 
tmpDist = zeros(1, n);
for i = 1 : 1 : size(Vec0, 1)
    tmpDist = tmpDist + Vec0(i, :).^2; 
end
tmpDist = sqrt(tmpDist); 
Vec1 = Vec0; 
for i = 1 : 1 : size(Vec0, 1)
    Vec1(i, :) = Vec1(i, :) ./ tmpDist; 
end
end

