function CvtPoseFun(DataRoot)
if nargin == 0
    DataRoot = 'J:\2017_11_26\Record-2017-11-26-08-35-45';
end
% DataRoot = 'G:\Data\XiAn\Record-2016-12-01-16-14-05(Est3RingR)\';
PoseDir = fullfile(DataRoot, 'Pose.txt');
PoseFid = fopen(PoseDir, 'r');
if PoseFid == -1
    fprintf('%s is not exist!', PoseDir);
    error('Stop!');
end
str = fullfile(DataRoot, 'SimplePose.txt');
fOutId = fopen(str, 'w');
nFrm = 0;
while 1
    tLine = fgetl( PoseFid );
    if tLine == -1
        break;
    end
    if isempty( tLine )
        continue;
    end
    [ wh wm ws wmm ] = strread( tLine, '---------%d:%d:%d:%d----------' );
    tLine = fgetl( PoseFid );
    if isempty(tLine)
        continue;
    end
    str = sprintf('%s %02d %02d %02d %03d\n', tLine, wh, wm, ws, wmm);
    fprintf(fOutId, str);
    if mod(nFrm, 100) == 0
        fprintf('Pose Converter, nFrm = %04d\n', nFrm);
    end
    nFrm = nFrm + 1;
end
fclose(fOutId);
fclose(PoseFid);

%%%%%%%%%%% Inserted Pose.
InsPoseDir = fullfile(DataRoot, 'InsertedPose.txt');
PoseFid = fopen(InsPoseDir, 'r');
if PoseFid == -1
    fprintf('%s is not exist!', InsPoseDir);
    error('Stop!');
end
str = fullfile(DataRoot, 'SimpleInsertedPose.txt');
fOutId = fopen(str, 'w');
nFrm = 0;
while 1
    tLine = fgetl( PoseFid );
    if tLine == -1
        break;
    end
    if isempty( tLine )
        continue;
    end
    for i = 1 : 1 : 3
        tLine = fgetl( PoseFid );
        str = sprintf('%s\n', tLine);
        fprintf(fOutId, str);
    end
    if mod(nFrm, 100) == 0
        fprintf('Inserted Pose Converter, nFrm = %04d\n', nFrm);
    end
    nFrm = nFrm + 1;
end
fclose(fOutId);
fclose(PoseFid);
end