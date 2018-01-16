function ScanData = CalScanData(IDX)
IdArray = unique(IDX);
ScanData = [];
for k = 1 : 1 : length(IdArray)
    Id = IdArray(k);
    if Id <= 0
        continue;
    end
    idx = find(IDX == Id);
    ScanData(:, end+1) = [Id; length(idx)];
end
end