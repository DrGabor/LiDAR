function InsPoseData = IterpPoseFun(PoseData, InsTime )
InsPoseData = [];
PoseTime = GetPoseFun(PoseData, 'time');
PoseTimeAll = GetTimeFun(PoseTime);
STimeAll = GetTimeFun(InsTime);
LPose = GetPoseFun(PoseData, 'local');
GPose = GetPoseFun(PoseData, 'global'); 
for id = 1 : 1 : length(STimeAll)
    if id == 7387
        bTest = 1; 
    end
    t = STimeAll(id);
    idx = find(t>=PoseTimeAll);
    idx0 = idx(end);
    idx = find(t<=PoseTimeAll);
    idx1 = idx(1);
    if idx0 == idx1
        s = 1.0; 
    else
        s = ( t - PoseTimeAll(idx0) ) / (PoseTimeAll(idx1) - PoseTimeAll(idx0) );
    end
    tmp = [];
    l = myInterp(LPose(idx0, :),  LPose(idx1, :), s);
    g = myInterp(GPose(idx0, :),  GPose(idx1, :), s); 
    InsPoseData(end+1, :) = [l g InsTime(id, :)];
    if mod(id, 100) == 0 
        str = sprintf('InterpPose, nFrm = %04d/%04d', id, length(STimeAll));
        disp(str); 
    end
    
end
end

%%% if angle difference is bigger than 180 degrees, the interpolated
%%% angle would point to a wrong direction, normalize the angle to
%%% avoid this situation.
function interpPose = myInterp(pose0, pose1, s)
tmp = pose0(4:6) - pose1(4:6);
NEffidx = find(abs(tmp) > 180.0);

for i = 1 : 1 : length(NEffidx)
    idx = NEffidx(i);
    Ang0 = pose0(3+i);
    Ang1 = pose1(3+i);
    if Ang0 < 0
        Ang1 = Ang1 - 360;
    end
    if Ang0 > 0
        Ang1 = Ang1 + 360;
    end
    pose1(3+i) = Ang1;
end
interpPose = (1-s)*pose0 + s*pose1;
end

