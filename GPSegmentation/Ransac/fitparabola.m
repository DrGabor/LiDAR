function [ V ] = fitparabola( XY )
    V = polyfit( XY(1, :), XY(2, :), 2 );
end

