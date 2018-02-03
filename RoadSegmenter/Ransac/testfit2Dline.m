% TESTFITLINE - demonstrates RANSAC line fitting
%
% Usage: testfit2Dline(outliers, sigma, t, feedback)
%
% Arguments:
%               outliers - Fraction specifying how many points are to be
%                          outliers.
%               sigma    - Standard deviation of inlying points from the
%                          true line.
%               t        - Distance threshold to be used by the RANSAC
%                          algorithm for deciding whether a point is an
%                          inlier. 
%               feedback - Optional flag 0 or 1 to turn on RANSAC feedback
%                          information.
%
%  Try using:  testfit2Dline(0.3, 0.05, 0.05)
%
% See also: RANSACFITPLANE, FITPLANE

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

% August 2006  testfitline created from testfitplane
%              author: Felix Duvallet

function testfit3Dline(outliers, sigma, t, feedback)

    close all;
    if nargin == 0
        outliers = 0.3;
        sigma = 0.05;
        t = 0.05;
        feedback = 1;
    end
    if nargin == 3
        feedback = 0;
    end
    
    % Hard wire some constants - vary these as you wish
    
    npts = 100;  % Number of 3D data points	
    
    % Define a line:
    %    Y = m*X
    %    Z = n*X + Y + b
    % This definition needs fixing, but it works for now
    
    m = 6;
    n = -3;
    b = -4;
    
    outsigma = 30*sigma;  % outlying points have a distribution that is
                          % 30 times as spread as the inlying points
    
    vpts = round((1-outliers)*npts);  % No of valid points
    opts = npts - vpts;               % No of outlying points
    
    % Generate npts points in the line
    X = rand(1,npts);
    
    Y = m*X + b;
    
    XY =  [X; Y];

    % Add uniform noise of +/-sigma
    XY = XY + (2*rand(size(XY))-1)*sigma;
    
    % Generate opts random outliers
    
    n = length(XY);
    ind = randperm(n);  % get a random set of point indices
    ind = ind(1:opts);  % ... of length opts
    
    % Add uniform noise of outsigma to the points chosen to be outliers.  
    XY(:,ind) = XY(:,ind)  +   sign(rand(2,opts)-.5).*(rand(2,opts)+1)*outsigma;    

    
    % Perform RANSAC fitting of the line
    [V, P, inliers] = ransacfit2Dline(XY, t, feedback);
    
    if(feedback)
        disp(['Number of Inliers: ' num2str(length(inliers)) ]);
    end

    % We want to plot the inlier points blue, with the outlier points in
    % red.  In order to do that, we must find the outliers.
    % Use setxor on all the points, and the inliers to find outliers
    %  (plotting all the points in red and then plotting over them in blue
    %  does not work well)
    oulier_points = setxor(transpose(XY), transpose(XY(:, inliers)), 'rows');
    oulier_points = oulier_points';
    
    % Display the cloud of outlier points
    figure(1); clf
    hold on;
    plot(oulier_points(1,:),oulier_points(2,:) , 'r*');

    % Plot the inliers as blue points
    plot(XY(1,inliers), XY(2, inliers), 'b*');

    % Display the line formed by the 2 points that gave the
    % line of maximum consensus as a green line
    line(P(1,:), P(2,:), 'Color', 'green', 'LineWidth', 4);
    
    %Display the line formed by the covariance fitting in magenta
    line(V(1,:), V(2, :), 'Color', 'magenta', 'LineWidth', 5);
    grid('on');
    % box('on'); 
    
    
