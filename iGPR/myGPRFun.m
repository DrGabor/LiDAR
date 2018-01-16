function [ymu, yms2] = myGPRFun(hyp, x, y, xs)
IS_SHOW = 0; 
if nargin == 0
    hyp = [];
    hyp.mean = [0.1 0.2];
    hyp.cov = log([1.2 1.5]);
    hyp.lik = log(0.1);
    x = [-5:0.1:5]';
    [Kff, ~] = CalKFun(x, x, hyp.cov);
    Mx = hyp.mean(1)*x + hyp.mean(2); 
    y = mvnrnd(Mx, Kff); 
    y = y'; 
    xs = [-10:0.1:10]'; % 10*rand(50, 1)-5;
    IS_SHOW = 1; 
end
m0 = hyp.mean(1);
m1 = hyp.mean(2);
sn = exp(hyp.lik);
sn2 = sn*sn;
n = length(x);
p = hyp.cov;
[Kuu, ~] = CalKFun(xs, xs, p);
[Kuf, ~] = CalKFun(xs, x, p);
[Kff, ~] = CalKFun(x, x, p);
A = Kff + sn2*eye(n); 
L = chol(A, 'lower');
iL = inv(L);
iA = iL'*iL;
% tt = diag(iA*A); 
% iA = inv(A); 
ymu = m0*xs + m1 + Kuf*iA*(y - m0*x-m1);
m = length(xs);
fs2 = diag(Kuu - Kuf*iA*Kuf');
yms2 = fs2+sn2;
bTest = 1;
if IS_SHOW
    Mu = ymu;
    S2 = yms2;
    Xs = xs;
    tData = 2.0; 
    a = Mu+tData*sqrt(S2); 
    b = flipdim(Mu-tData*sqrt(S2),1); 
    F = [a;b];
    figure;
    hold on;
    grid on;
    axis equal;
    %%%%%%%%%% draw envelop.
%     plot(xs, a, 'r.--'); 
%     plot(xs, b, 'g.--'); 
    fill([Xs; flipdim(Xs,1)], F, [7 7 7]/8)
    plot(Xs, Mu, 'b.');
    plot(x, y, 'r.'); 
end
end
