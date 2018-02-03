function [Coor, EffIdx, ind] = CvtPtsToVoxelFun( data, minVal, Size, gridRes )
Coor = [];
EffIdx = [];
ind = [];
Dim = size(data, 1); 
for i = 1 : 1 : Dim
    Coor(i, :) = ceil(( data(i, :) - minVal(i) )/gridRes);
end
if Dim == 3
    EffIdx = find( Coor(1, :) >= 1 & Coor(1, :) <= Size(1) & ...
        Coor(2, :) >= 1 & Coor(2, :) <= Size(2) & ...
        Coor(3, :) >= 1 & Coor(3, :) <= Size(3) );
    ind = sub2ind(Size, Coor(1, EffIdx), Coor(2, EffIdx), Coor(3, EffIdx));
end
if Dim == 2
    EffIdx = find( Coor(1, :) >= 1 & Coor(1, :) <= Size(1) & ...
        Coor(2, :) >= 1 & Coor(2, :) <= Size(2) );
    ind = sub2ind(Size, Coor(1, EffIdx), Coor(2, EffIdx));
end
end