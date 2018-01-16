function [ varargout ] = TrimmedICP(Md, MovData, RefData, Tf0, MaxIter,TrMin,TrMax,lamda)
if nargin == 0
    clc; close all;
    %load ../../StanfordData/bunny;
    RefData = rand(3, 100); %bunny{1}';
    MovData = rand(3, 100);% bunny{2}';
    Md = createns(RefData');
    Tf0 = [eye(3) zeros(3, 1)]; 
    MaxIter = 50; 
    TrMin = 0.35; 
    TrMax = 1.0; 
    lamda = 2.0; 
end
Dim = size(MovData, 1); 
PreMSE= 10^5;   CurMSE= 10^6;  Iter = 1; 
RelErr = 1.0;
R0 = Tf0(1:Dim, 1:Dim); 
T0 = Tf0(1:Dim, end); 
Tf = [R0 T0; zeros(1, Dim) 1 ]; 
for Iter = 1 : 1 : MaxIter       % (abs(CurMSE-PreMSE)>10^(-9))   % this threshold is sensitive to resolutions....for bunny, this value should be 1e-12.
    TData = Loc2Glo(MovData, Tf(1:Dim, 1:Dim)', Tf(1:Dim, end) );
    [corr,TD] = knnsearch( Md,TData(1:Dim, :)');
    SortTD2 = sortrows(TD.^2); % Sort the correspongding points
    minTDIndex = floor(TrMin*length(TD)); % Get minimum index of TD
    maxTDIndex = ceil(TrMax*length(TD)); % Get maxmum index of TD    
    TDIndex = [minTDIndex : maxTDIndex]';
    mTr = TDIndex./length(TD);
    mCumTD2 = cumsum(SortTD2);
    mMSE = mCumTD2(minTDIndex : maxTDIndex)./TDIndex;
    mPhi = ObjectiveFunction(mMSE, mTr);  
    PreMSE=CurMSE;
    [CurMSE, nIndex] = min(mPhi);    
    Trim = mTr(nIndex); % Update Tr for next step    
    corr(:,2) = [1 : length(corr)]';
    % Sort the corresponding points
    corrTD = [corr, TD];
    SortCorrTD = sortrows(corrTD, 3);
    
    TrLength = floor(Trim*size(SortCorrTD,1)); % The number of corresponding points after trimming
    TCorr = SortCorrTD(1:TrLength, 1:2);     % Trim the corresponding points according to overlap parameter Tr
    % Register MData with TData
    dM = RegFun(RefData(:, TCorr(:, 1)), TData(:, TCorr(:, 2)) ); 
    % dM = reg(RefData, TData, TCorr);
    Tf = dM * Tf; 
    % [M, TCorr, scan] = CalRtPhi(model, data, SortCorrTD, Trim);
    Err = [norm(dM(1:Dim, 1:Dim)-eye(Dim)) norm(dM(1:Dim, end))];
    if max(Err) <= 1e-5
        break;
    end
    RelErr = (PreMSE - CurMSE) / abs(PreMSE); 
    if abs(RelErr) <= 1e-3
        break; 
    end
end
IS_SHOW = 0;
if IS_SHOW
    Res = CalRes(MovData); 
    h = figure; 
    ICP_PlotFun(MovData, RefData, Tf, Res, h); 
    title('Trimmed ICP');
end;
if nargout == 1
    varargout{1} = Tf; 
end
if nargout == 2 
    varargout{1} = Tf(1:Dim, 1:Dim); 
    varargout{2} = Tf(1:Dim, end); 
end
bTest = 1; 

end

%%%%%%%%%%%%%%%%%%%%Integrated Function%%%%%%%%%%%%%%%%%%%%
% %% Calculate R,t,Phi based on current overlap parameter
% function [M,TCorr,TData] = CalRtPhi(Model, scan, SortCorrTD,Tr)
% 
% TrLength = floor(Tr*size(SortCorrTD,1)); % The number of corresponding points after trimming
% TCorr = SortCorrTD(1:TrLength, 1:2);     % Trim the corresponding points according to overlap parameter Tr
% % Register MData with TData
% [M] = reg(Model(1:3,:), scan(1:3,:), TCorr);
% % To obtain the transformation data
% TData = M*scan;
% 
% end
% 
% function [M,TCorr,TData] = CalRtPhi(Model, scan, SortCorrTD,Tr)
% 
% TrLength = floor(Tr*size(SortCorrTD,1)); % The number of corresponding points after trimming
% TCorr = SortCorrTD(1:TrLength, 1:2);     % Trim the corresponding points according to overlap parameter Tr
% % Register MData with TData
% [M] = reg(Model(1:3,:), scan(1:3,:), TCorr);
% % To obtain the transformation data
% TData = M*scan;
% 
% end
%%%%%%%%%%%%%%% Calculate the registration matrix %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% T(TData)->MData %%%%%%%%%%%%%%%%%%%%%%%%%
% SVD solution
function [M] = reg(Model, Data, corr)

n = length(corr); 
M = Model(:,corr(:,1)); 
mm = mean(M,2);
S = Data(:,corr(:,2));
ms = mean(S,2); 
Sshifted = [S(1,:)-ms(1); S(2,:)-ms(2); S(3,:)-ms(3)];
Mshifted = [M(1,:)-mm(1); M(2,:)-mm(2); M(3,:)-mm(3)];
K = Sshifted*Mshifted';
K = K/n;
[U A V] = svd(K);
R1 = V*U';
if det(R1)<0
    B = eye(3);
    B(3,3) = det(V*U');
    R1 = V*B*U';
end
t1 = mm - R1*ms;
M=[];
M(1:3,1:3)=R1;
M(1:3,4)=t1;
M(4,:)=[0,0,0,1];
end

function [Phi] = ObjectiveFunction(MSE, TrB)
lamga= 2;
Phi = MSE./((TrB).^((1+lamga)));
end


