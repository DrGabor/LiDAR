function showCorresFunc(meshX,meshY, corIdx1,corIdx2, offset)

N = length(corIdx1);
cmap = hsv(max(N,1)); 
meshX.vertices = meshX.vertices - repmat(offset,length(meshX.vertices),1);
figure; trisurf(meshX.faces,meshX.vertices(:,1),meshX.vertices(:,2),meshX.vertices(:,3)); axis equal; axis image; shading interp; lighting phong; view([0,0]);  camlight left; hold on; colormap([1,1,0]*0.9);
hold on;
trisurf(meshY.faces,meshY.vertices(:,1),meshY.vertices(:,2),meshY.vertices(:,3)); axis equal; axis image; shading interp; lighting phong; view([0,0]);   hold on; colormap([1,1,1]*0.9);
% show correspondence                        
for m = 1:N       
        line([meshX.vertices(corIdx1(m),1), meshY.vertices(corIdx2(m),1)]',...
               [meshX.vertices(corIdx1(m),2), meshY.vertices(corIdx2(m),2)]',...
               [meshX.vertices(corIdx1(m),3), meshY.vertices(corIdx2(m),3)]','Color',cmap(m,:),'LineWidth',0.5);
end
axis image, axis off,


