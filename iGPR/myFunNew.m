function [f, df] = myFunNew(X0, x, y, sn)
if nargin == 0
    clc; close all;
    x = [-20:0.1:60]';
    y = rand(length(x), 1);
    n = length(x);
    sn = 0.1;
    X0 = [0.1 0.1 log([2.0 2.0])]';
end
if isrow(x)
    x = x'; 
end
if isrow(y)
    y = y'; 
end
a = X0(1);
b = X0(2);
l = exp(X0(3));
l2 = l*l;
sf = exp(X0(4));
sf2 = sf*sf;
sn2 = sn^2; 

n = length(x); 
[K0, dK] = CalKFun(x, x, X0(3:4)); 
K = K0 + sn2*eye(n); 
scale = 1000; 
tt = eig(scale*K)/scale;
%%%%%%%%% make sure K is positive definite.
if min(tt) < 0.0
    K = K - min(tt)*eye(n); 
end
if rank(K) < n
    K = K + 1e-6*eye(n); 
end
% save('./K.mat'); 
L = chol(K, 'lower');
iL = inv(L); 
iK = iL'*iL; 
logDetK = 2*sum(log(diag(L))); 

mx = a*x + b; 
f = -0.5*(y-mx)'*iK*(y-mx) - 0.5*logDetK-0.5*n*log(2*pi);
%%%%%%%%%% calculate gradient. 

Alpha = iK*(y-mx);
df2m = [x'; ones(1, n)] * Alpha;
A = Alpha * Alpha' - iK;
df2c = zeros(2, 1); 
for i = 1 : 1 : 2 
     nSum = 0; 
     for id = 1 : 1 : n
         nSum = nSum + A(id, :) * dK{i}(:, id); 
     end
     df2c(i) = 0.5*nSum; 
end
df = [df2m; df2c]; 
f = -f; 
df = -df; 
bTest = 1; 
end