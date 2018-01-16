function [ AngleNew ] = wrapTo2PiFun(Angle)
AngleNew = Angle; 
idx = find(Angle < 0 );
AngleNew(idx) = AngleNew(idx) + 2*pi; 
end

