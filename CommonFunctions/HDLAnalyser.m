%% convert smaller binary file .txt(~298kB) to 3D point cloud, this function costs about 60~90ms.
function pcData = HDLAnalyser( DataDir, R, T )
if nargin == 0 
    DataDir = 'F:\Data\Record-2016-10-24-09-49-36(RealTrafficScene)\BinaryData\Binary100.txt'; 
    R = eye(3);
    T = [1.2; 0.0; 0.0];
end
if nargin == 1
    R = eye(3);
    T = [1.2; 0.0; 0.0];
end
%% Correction parameter rectified.
tmp = load('C:\Program Files\MATLAB\R2015a\bin\Velodyne\CorrPara.txt'); % LaserId  RotCorrAng  VertCorrAng  DistCorr  VertOffset HorizOffset;
SeqMap = [36 37 58 59  38 39 32 33 40 41 34 35 48 49 42 43 50 51 44 45 52 53 46 47 60 61 54 55 62 63 56 57 4  5  26 27 6  7  0  1  8 9 2  3  16 17 10 11 18 19 12 13 20 21 14 15 28 29 22 23 30 31 24 25];
[~, Idx] = sort(SeqMap);
tmp = tmp( Idx,: );
Para = tmp(:, 2:end );   % RotCorrAng  VertCorrAng DistCorr VertOffset HorizOffset;
CosSinTab = [cosd( (0:1:3999) * 0.09 );
    sind( (0:1:3999) * 0.09) ];
%% Read binary data from .txt.
fid = fopen(DataDir, 'rb' );
BinaryData = fread(fid, 'uint8=>uint16');
UDPNum = length(BinaryData) / 1206;
A = reshape( BinaryData, 1206, UDPNum);
fclose(fid);
%% Extract raw angle, distance and intensity.
RawDistData = zeros(3, UDPNum * 6, 64);
for id = 1 : 1 : 6
    Idx = (1 + (id-1)*200) : 1 : id * 200;
    Buffer = A(Idx, :);
    
    UpLow1 = Buffer(1, :) + 256 * Buffer(2, :);
    UpLow2 = Buffer(101, :) + 256 * Buffer(102, :);
    Idx = find(UpLow1 ~= hex2dec('eeff') & UpLow2 ~= hex2dec('ddff') );
    if ~isempty(Idx)
        disp('The Up block(0xEEFF) and the Low Block(0xDDFF) is not right!\n');
    end
    RotAng1 = Buffer(3, :) + 256 * Buffer(4, :);
    RotAng2 = Buffer(103, :) + 256 * Buffer(104, :);
    Idx = find( (RotAng1 == RotAng2) ~= 1 );
    if ~isempty(Idx)
        disp('Two Angle in the Block is not equal! \n');
    end
    AngID = floor(RotAng1 / 9);
    Idx = ( 1 + (id-1) * UDPNum ) : 1 : id * UDPNum;
    for i = 0 : 1 : 31
        Id0 = SeqMap(i+1)+1;
        Id1 = SeqMap(i+33)+1;
        RawDistData(1, Idx, Id0)  = AngID;
        RawDistData(2, Idx, Id0)  = Buffer(i*3+5, :) + 256 * Buffer(i*3+6, :);
        RawDistData(3, Idx, Id0)  = Buffer(i*3+7, :);
        RawDistData(1, Idx, Id1) = AngID;
        RawDistData(2, Idx, Id1) = Buffer(i*3+105, :) + 256 * Buffer(i*3+106, :);
        RawDistData(3, Idx, Id1) = Buffer(i*3+107, :);
    end
end
%% median filter. it seems this will introduce some noise.....
USE_MEDIAN_FILTER = 0;
if USE_MEDIAN_FILTER
   for i = 1 : 1 : 64
    Dist = RawDistData(2, :, i);
    % RawDistData(2, :, i) = medfilt1( Dist, n );
    tmp = [ Dist( (end-n+1):1:end) Dist Dist(1:1:n) ];
    tmp = medfilt1( tmp, n );
    RawDistData(2, :, i) = tmp((n+1):1:end-3);
end 
end
n = 3;

%% analysis data.
pcData = [];
Radius = [];
RotAngID = [];
for i = 1 : 1 : 64
    RotAngID = RawDistData(1, :, i);
    Radius   = RawDistData(2, :, i);
    Idx0 = find(Radius > 0.0 );
    RotAngID = RotAngID(Idx0);
    Radius   = Radius(Idx0);
    Gama = Para(i, 1);  % rotational angle.
    Beta = Para(i, 2);  % vertical angle.
    D = Para(i, 3);
    V = Para(i, 4) / 100.0;
    H = Para(i, 5) / 100.0;
    dDist = ( 0.2 * Radius + D ) / 100.0;
    RotAng = RotAngID * 0.09 - Gama;   % 0.09 is the angle's resolution.
    
    c = CosSinTab(1, RotAngID+1 );
    s = CosSinTab(2, RotAngID+1 );
    cRotAng = c * cosd(Gama) + s * sind(Gama);
    sRotAng = s * cosd(Gama) - c * sind(Gama);
    dDistXY = dDist * cosd(Beta) - V * sind(Beta);
    x =  dDistXY .* cRotAng + H * sRotAng;
    y = -dDistXY .* sRotAng + H * cRotAng;
    z = dDist * sind(Beta) + V * cosd(Beta);
    tmp = [x; y; z; i * ones(1, length(x))];
    Idx = find( dDistXY >= 2.0 & dDistXY <= 100.0 );
    pcData = [pcData tmp(:, Idx)];
end
pcData(1:3, :) = R * pcData(1:3, : ) + repmat(T, 1, size(pcData, 2) );
end

