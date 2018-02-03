function varargout = VoxelizeFun( data, gridRes )
out = [];
ReduceData = [];
CoorI = [];
Res = 1e-6;
minVal = min(data') - Res;
maxVal = max(data') + Res;
Size = ceil((maxVal - minVal)/gridRes);
Coor_Raw = [];
Dim = size(data, 1); 
for i = 1 : 1 : Dim
    Coor_Raw(i, :) = ceil(( data(i, :) - minVal(i) )/gridRes);
end
if Dim == 3
    ind = sub2ind(Size, Coor_Raw(1, :), Coor_Raw(2, :), Coor_Raw(3, :));
end
if Dim == 2
    ind = sub2ind(Size, Coor_Raw(1, :), Coor_Raw(2, :));
end
[B, I] = sort(ind);
[C, ia, ic] = unique(B);
if Dim == 3
    [X Y Z] = ind2sub(Size, C);
    CoorI = [X; Y; Z];
end
if Dim == 2
    [X Y] = ind2sub(Size, C);
    CoorI = [X; Y];
end
parfor i = 1 : 1 : length(ia)
    if i == 1
        bTest = 1;
    end
    id0 = ia(i);
    if i == length(ia)
        id1 = length(ind);
    else
        id1 = ia(i+1)-1;
    end
    idx = I(id0:id1);
    tmp = ind(idx);
    a = unique(tmp);
    if length( unique(tmp) ) ~= 1 | a ~= C(i)
        error('Error!');
    end
    if length(idx) > 1
        pt = mean( data(:, idx)')';
    else
        pt = data(:, idx);
    end
    ReduceData = [ReduceData pt];
    bTest = 1;
end
GridData = CoorI * gridRes;
GridData = bsxfun(@plus, GridData, minVal'-gridRes/2*ones(Dim, 1));

if nargout == 1
    out.reduceData = ReduceData;
    out.coorI = CoorI;
    out.gridData = GridData;
    varargout{1} = out;
end
if nargout == 3 
    varargout{1} = ReduceData; 
    varargout{2} = CoorI; 
    varargout{3} = GridData;
end
% [X Y Z] = ind2sub(Size, C);
% tmp = [X; Y; Z] * gridSize;
% tmp = bsxfun(@minus, tmp, gridSize/2);
% %%%%%%%%%%%% ReduceData is more accurate than ReduceGridData.
% ReduceGridData = bsxfun(@plus, tmp, minVal');
end

