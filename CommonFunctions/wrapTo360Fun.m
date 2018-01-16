function [ AngleNew ] = wrapTo360Fun( Angle )
AngleNew = Angle; 
idx = find(Angle < 0 );
AngleNew(idx) = AngleNew(idx) + 360; 
end

