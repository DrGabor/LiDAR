function [K, dK] = CalKFun(x, y, p)
l = exp(p(1));
sf = exp(p(2));
l2 = l*l;
sf2 = sf*sf;
m = length(x);
n = length(y);
dK = {};
X = repmat(x, 1, n); 
Y = repmat(y', m, 1); 
R = (X - Y) .*(X - Y);  
K = sf2*exp(-0.5/l2*R); 
dK{1} = 1/l2 * R.*K;
dK{2} = 2*K;
end
