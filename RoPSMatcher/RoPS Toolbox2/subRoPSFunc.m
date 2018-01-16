function subRoPS =  subRoPSFunc(projNeighbor,binSize)
%get the distribution matrix
neighbNum = length(projNeighbor);
distrMatrix = zeros(binSize,binSize);
minX = min(projNeighbor(:,1));
stepX = (max(projNeighbor(:,1))-minX)/binSize;
minY = min(projNeighbor(:,2));
stepY = (max(projNeighbor(:,2))-minY)/binSize;

if stepX==0 || stepY==0
    subRoPS = [0,0,0,0,0];
    return;
end

for k=1:neighbNum
    idxX = ceil((projNeighbor(k,1) - minX)/stepX);
    idxY = ceil((projNeighbor(k,2) - minY)/stepY);
    if idxX>binSize      idxX = binSize;  end
    if idxX<1                idxX = 1;            end
    if idxY>binSize     idxY = binSize;   end
    if idxY<1                idxY = 1;            end
    distrMatrix(idxX,idxY) = distrMatrix(idxX,idxY)+1;
end
distrMatrix = distrMatrix/neighbNum;%normalization
%calculate the moment of this distribution matrix
meanX = 0;
meanY = 0;
pde = 0;
for idxX = 1:binSize
    for idxY = 1:binSize
        meanX = meanX+idxX*distrMatrix(idxX,idxY);
        meanY = meanY+idxY*distrMatrix(idxX,idxY);
        if distrMatrix(idxX,idxY)>0
            pde = pde - distrMatrix(idxX,idxY)*log2(distrMatrix(idxX,idxY));
        end
    end
end
u11 = 0;
u21 = 0;
u12 = 0;
u22 = 0;

for idxX = 1:binSize
    for idxY = 1:binSize
        u11 = u11+(idxX-meanX)*(idxY-meanY)*distrMatrix(idxX,idxY);
        u21 = u21+(idxX-meanX)^2*(idxY-meanY)*distrMatrix(idxX,idxY);
        u12 = u12+(idxX-meanX)*(idxY-meanY)^2*distrMatrix(idxX,idxY);
        u22 = u22+(idxX-meanX)^2*(idxY-meanY)^2*distrMatrix(idxX,idxY);
    end
end 
 subRoPS = [u11,u21,u12,u22,pde];