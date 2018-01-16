function [TimeAll] = GetTimeFun(Time)
if size(Time, 2) ~= 4
    error('Time error!\n'); 
end
TimeAll = ( ( Time(:, 1)*60 + Time(:, 2) ) * 60 + Time(:, 3) )*1000 + Time(:, 4); 
end

