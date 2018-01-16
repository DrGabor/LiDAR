function out = LaserEmitFun(data, ExLen, Res0, IS_SHOW)
if nargin == 0
    clc; close all;
    Tao0 = deg2rad(10.0);   %%%%%%%%%%% initial tangible angle.
    Kapa = 1/650;     %%%%%%%%%%% initial curvature. According to highway standard, for moutainous region, R >= 250m; for plain, R >=650m.
    c = 4*1e-4; 0.001;   %%%%%%%%%% change of curvature.
    f_cos = @(x, c0, c1, c2) cos(c0 + c1*x+c2*x.^2/2);
    f_sin = @(x, c0, c1, c2) sin(c0 + c1*x+c2*x.^2/2);
    x0 = 0.0;
    y0 = 0.0;
    K_g = [];
    Ang_g = [];
    x_g = [];
    y_g = [];
    for l = 0.0:0.1:20
        x = x0 + integral(@(x) f_cos(x, Tao0, Kapa, c), 0, l );
        y = y0 + integral(@(x) f_sin(x, Tao0, Kapa, c), 0, l );
        K_g(end+1) = Kapa + c*l;
        Ang_g(end+1) = Tao0 + Kapa*l + c*l^2/2;
        x_g(end+1) = x;
        y_g(end+1) = y;
    end
    x_obv = x_g;
    y_obv = y_g; %  + 0.1*( rand(1, length(y_g)) - 0.5 );
    data = [x_obv; y_obv];
    ExLen = 10;
    Res0 = 0.5;
    IS_SHOW = 1;
end
x_obv = data(1, :);
y_obv = data(2, :);
if IS_SHOW
    figure;
    hold on;
    grid on;
    plot(x_obv, y_obv, 'b.');
end

%%%%%% for data, the first and last point is out of consideration.
PtsVec = data(:, 2:end-1);
FrontVec = data(:, 3:end) - data(:, 2:end-1);
BackkVec = data(:, 1:end-2) - data(:, 2:end-1);
FrontVec = NormalizeVecFun(FrontVec);
BackkVec = NormalizeVecFun(BackkVec);
VecDiff = FrontVec - BackkVec;
NormDiff = NormalizeVecFun(VecDiff);
Ang = 90.0;
RL = [cosd(Ang) -sind(Ang)
    sind(Ang)  cosd(Ang) ];
Ang = -90.0;
RR = [cosd(Ang) -sind(Ang)
    sind(Ang)  cosd(Ang) ];
ss = struct('L', [], 'R', [] );
out = [];
Range = Res0:ExLen;
RangeArray = [Range; Range];
nLen = length(Range);
for i = 1 : 1 : size(NormDiff, 2)
    tmp = ss;
    Pts = RangeArray .* repmat( NormDiff(:, i), 1, nLen );
    pt0 = data(:, i+1);
    tmp.L = Loc2Glo(Pts, RL', pt0);
    tmp.R = Loc2Glo(Pts, RR', pt0);
    out = [out tmp];
end

% VecDiff = tmp * NormalizeVecFun(VecDiff);

LVec = Loc2Glo(VecDiff, RL', [0; 0]);
RVec = Loc2Glo(VecDiff, RR', [0; 0]);
if IS_SHOW
    figure;
    hold on;
    grid on;
    axis equal;
    plot(PtsVec(1, :), PtsVec(2, :), 'k.');
    % quiver(PtsVec(1, :), PtsVec(2, :), BackkVec(1, :), BackkVec(2, :), 0.05, 'Color', 'b');
    quiver(PtsVec(1, :), PtsVec(2, :), LVec(1, :), LVec(2, :), 0.05, 'Color', 'b');
    quiver(PtsVec(1, :), PtsVec(2, :), RVec(1, :), RVec(2, :), 0.05, 'Color', 'r');
    for i = 1 : 1 : length(out)
        L = out(i).L;
        R = out(i).R;
        plot(L(1, :), L(2, :), 'g.-');
        plot(R(1, :), R(2, :), 'r.-');
    end
end
bTest = 1;
end
