% Usage  [L, inliers] = ransacfitline(XYZ, t, feedback)
%
% Arguments:
%          XYZ - 3xNpts array of xyz coordinates to fit line to.
%          t   - The distance threshold between data point and the line
%                used to decide whether a point is an inlier or not.
%          feedback - Optional flag 0 or 1 to turn on RANSAC feedback
%                     information.
%
% Returns:.
%           V - Line obtained by a simple fitting on the points that
%               are considered inliers.  The line goes through the
%               calculated mean of the inlier points, and is parallel to
%               the principal eigenvector.  The line is scaled by the
%               square root of the largest eigenvalue.
%               This line is a n*2 matrix.  The first column is the
%               beginning point, the second column is the end point of the
%               line.
%           L - The two points in the data set that were found to
%               define a line having the most number of inliers.
%               The two columns of L defining the two points.
%           inliers - The indices of the points that were considered
%                     inliers to the fitted line.
%
% See also:  RANSAC, FITPLANE, RANSACFITPLANE

% Copyright (c) 2003-2006 Peter Kovesi and Felix Duvallet (CMU)
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Aug  2006 - created ransacfitline from ransacfitplane
%             author: Felix Duvallet

function [V, L, inliers] = ransacfitparabola(XY, t, feedback)
    
    if nargin == 2
	feedback = 0;
    end
    
    [rows, npts] = size(XY);
    
    if rows ~= 2
        error('data is not 2D');
    end
    
    if npts < 2
        error('too few points to fit line');
    end
    
    s = 3;  % Minimum No of points needed to fit a line.
        
    fittingfn = @defineParabola;
    distfn    = @paraptdist;
    degenfn   = @isdegenerate;

    [L, inliers] = ransac(XY, fittingfn, distfn, degenfn, s, t, feedback);
    
    V = fitparabola(XY(:, inliers));
    
%------------------------------------------------------------------------
% Function to define a parabola given 3 data points as required by
% RANSAC.

function L = defineParabola(X)
    L = X;
%------------------------------------------------------------------------
% Function to calculate distances between a line and an array of points.
function [inliers, L] = paraptdist(L, X, t)
    p = polyfit( L(1, :), L(2, :), 2 );
    YEst = polyval(p, X(1, :) );
    d = YEst - X(2, :);
    inliers = find(abs(d) < t);
    
%------------------------------------------------------------------------
% Function to determine whether a set of 2 points are in a degenerate
% configuration for fitting a line as required by RANSAC.
% In this case two points are degenerate if they are the same point
% or if they are exceedingly close together.

function r = isdegenerate(X)
    %find the norm of the difference of the two points
    % this will be 0 iff the two points are the same (the norm of their
    % difference is zero)
    r0 = norm(X(:,1) - X(:,2)) < eps;
    r1 = norm(X(:,1) - X(:,3)) < eps;
    r = r0 | r1;