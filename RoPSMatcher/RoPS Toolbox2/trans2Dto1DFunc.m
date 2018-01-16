function oneDfeature = trans2Dto1DFunc(twoDfeature)
rowSize = size(twoDfeature,1);
oneDfeature = [];
for i=1:rowSize
    oneDfeature = [oneDfeature,twoDfeature(i,:)];
end
