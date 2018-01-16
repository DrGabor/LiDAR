function [f, df] = mySparseFun(X0, x, y, xu, s, sn)
if nargin == 0
    clc; close all;
    x = [-20:0.1:60]';
    y = rand(length(x), 1);
    xu = [-20:0.1:60]';
    n = length(x);
    s = 1.0;
    sn = 0.1;
    X0 = [0.1 0.1 log([2.0 2.0])]';
end
a = X0(1);
b = X0(2);
l = exp(X0(3));
l2 = l*l;
sf = exp(X0(4));
sf2 = sf*sf;
sn2 = sn*sn; 

m = length(xu);
n = length(x); 
[Kuu, dKuu] = CalKFun(xu, xu, X0(3:4));
Kuu = Kuu + sf2*1e-2*eye(m);    % stablize. 
save('./KInfo.mat', 'Kuu', 'X0', 'xu'); 
if rank(Kuu) < m
    Kuu = Kuu + 1e-6*eye(m); 
end
L0 = chol(Kuu, 'lower'); 
iL0 = inv(L0); 
iKuu = iL0'*iL0;
logDetKuu = 2*sum(log(diag(L0))); 

[Kfu, dKfu] = CalKFun(x, xu, X0(3:4)); 
Kuf = Kfu'; 
dKuf = {dKfu{1}', dKfu{2}'}; 

Qff = Kfu*iKuu*Kuf;
DiagKff = sf2*eye(n); 
Gama = sn2*eye(n) + s*(DiagKff - diag(diag(Qff)));
iGama = diag( 1.0 ./ diag(Gama) ); 
logDetGama = sum( log(diag(Gama)) ); 
ApxK = Gama + Qff;   % use Woodbury to calculate its inverse and determinants: Z = Gama, U = V = Kfu, W = iKuu. 
S = Kuu + Kfu'*iGama*Kfu; 
% S = S + 1e-10*eye(m);  % for numerical stabitity. 
% if rank(S) < m
%     S = S + 1e-6*eye(m); 
% end
L1 = chol(S, 'lower'); 
logDetS = 2*sum( log(diag(L1)) ); 
iL1 = inv(L1); 
iS = iL1'*iL1; 
iApxK = iGama - iGama*Kfu*iS*Kfu'*iGama;
tmp = iApxK*ApxK; 
tt = diag(tmp); 
% max(tt) - min(tt)

logDetK = logDetGama - logDetKuu + logDetS;   % use Woodbury equation.
mx = a*x + b; 
f = -0.5*(y-mx)'*iApxK*(y-mx) - 0.5*logDetK-0.5*n*log(2*pi);
% ff = -0.5*(y-mx)'*pinv(ApxK)*(y-mx) - 0.5*log(det(ApxK)) - 0.5*n*log(2*pi); 
% -0.5*(y-mx)'*pinv(ApxK)*(y-mx)+0.5*(y-mx)'*iApxK*(y-mx)
%%%%%%%%%% calculate gradient. 
dQ = {}; 
dGama = {}; 
dDiagK = {zeros(n, n), 2*DiagKff}; 
dK = {}; 
Alpha = iApxK*(y-mx);
df2m = [x'; ones(1, n)] * Alpha;
A = Alpha * Alpha' - iApxK;
df2c = zeros(2, 1); 
for i = 1 : 1 : 2
     dQ{i} = dKfu{i}*iKuu*Kuf + Kfu*iKuu*dKuf{i}+Kfu*(-iKuu*dKuu{i}*iKuu)*Kuf;
     dGama{i} = s*dDiagK{i} - s*diag(diag(dQ{i})); 
     dK{i} = dQ{i} + dGama{i}; 
     nSum = 0; 
     for id = 1 : 1 : n
         nSum = nSum + A(id, :) * dK{i}(:, id); 
     end
     df2c(i) = 0.5*nSum; 
end
% TestWoodFun(Gama, iKuu, Kfu, Kfu); 
%%%%%%%%% Test dQ. 
% Eps = 1e-2; 
% for i = 1 : 1 : 2
%     v = zeros(2, 1); 
%     v(i) = Eps; 
%     XNew = X0(3:4) + v;
%     [KfuNew, ~] = CalKFun(x, xu, XNew);
%     [KuuNew, ~] = CalKFun(xu, xu, XNew); 
%     QffNew = KfuNew*inv(KuuNew)*KfuNew'; 
%     dQ_Est = (QffNew - Qff) / Eps; 
%     dKuuNew = (KuuNew - Kuu)/Eps; 
%     tmpDiff = dQ{i} - dQ_Est; 
%     norm(tmpDiff(:)) 
% end

df = [df2m; df2c]; 

f = -f; 
df = -df; 
bTest = 1; 
end

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
% tmp = zeros(m,n); 
% for i = 1 : 1 : m
%     for j = 1 : 1 : n
%         r = x(i) - y(j); 
%         r2 = r*r; 
%         tmp(i, j) = sf2*exp(-0.5*r2/l2); 
%     end
% end
% bTest = 1; 
end

function [K] = TestWoodFun(Z, W, U, V)
K = Z + U * W * V'; 
InvZ = pinv(Z); 
InvW = pinv(W); 
S = InvW+V'*InvZ*U; 
InvK = InvZ - InvZ*U*pinv(S)*V'*InvZ;
detK = det(Z)*det(W)*det(InvW + V'*InvZ*U); 

A = K*InvK; 
figure; 
tt = diag(A); 
plot(tt, 'b.' ); 
max(tt) - min(tt)
bTest = 1; 
end
