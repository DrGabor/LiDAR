%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function PlotClusterinResult(X, IDX)
    nDim = size(X, 2); 
    k=max(IDX);

    Colors=hsv(k);

    Legends = {};
    for i=0:k
        Xi=X(IDX==i,:);
        if i~=0
            Style = 'x';
            MarkerSize = 8;
            Color = Colors(i,:);
            Legends{end+1} = ['Cluster #' num2str(i)];
        else
            Style = 'o';
            MarkerSize = 6;
            Color = [0 0 0];
            if ~isempty(Xi)
                Legends{end+1} = 'Noise';
            end
        end
        if ~isempty(Xi)
            str = sprintf('Id=%02d, Len=%02d', i-1, length(Xi(:, 1))); 
            pt = Xi(1, :); 
            if nDim == 2 
                text(pt(1), pt(2), str, 'FontSize', 12 ); 
                plot(Xi(:,1),Xi(:,2),Style,'MarkerSize',MarkerSize,'Color',Color);
            else
                text(pt(1), pt(2), pt(3), str, 'FontSize', 12 );
                plot3(Xi(:,1),Xi(:,2),Xi(:, 3), Style,'MarkerSize',MarkerSize,'Color',Color);
            end
        end
        hold on;
    end
    hold off;
    axis equal;
    grid on;
    % legend(Legends);
    % legend('Location', 'NorthEastOutside');

end