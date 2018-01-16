function Data = GetPoseFun(RawPoseData, str)
%%%%% each row's format
% px py pz head pitch roll Speed gpsPx gpsPy gpsPz gpsHead gpsPitch gpsRoll gpsSpeed wh wm ws wmm
Data = []; 
if strcmp(str, 'local')
    Data = RawPoseData(:, 1:7); 
end
if strcmp(str, 'global')
    Data = RawPoseData(:, 8:14); 
end
if strcmp(str, 'time')
    Data = RawPoseData(:, end-3:end); 
end
end