clc; close all;
% x = gallery('uniformdata',[10,1],0);
% y = gallery('uniformdata',[10,1],1);
% DT = delaunayTriangulation(x,y);
% 
% k = convexHull(DT)
% figure
% plot(DT.Points(:,1),DT.Points(:,2), '.','markersize',10);
% hold on
% plot(DT.Points(k,1),DT.Points(k,2),'r')
% hold off

figure; 
x = dTreePtsL(1, :)';
y = dTreePtsL(2, :)';
% plot(x,y,'.')
DT = delaunayTriangulation([x y] );
triplot(DT);
axis equal;
% xlim([-0.2 1.2])
% ylim([-0.2 1.2])

% k = boundary(x,y);
% hold on;
% plot(x(k),y(k));