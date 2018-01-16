function [tmp] = PointsNorm( Pts )
    Dim = size(Pts, 1);
    tmp = [];
    if Dim == 2 
        tmp = Pts(1, :).^2 + Pts(2, :).^2;
    else
        tmp = Pts(1, :).^2 + Pts(2, :).^2 + Pts(3, :).^2;
    end
end

