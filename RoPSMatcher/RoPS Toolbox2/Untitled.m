clc; close all; 
CorrIdx0 = mesh1.keypntIdx(corPntIdx(:,1));
CorrIdx1 = mesh2.keypntIdx(corPntIdx(:,1));
RefData = mesh1.vertices( CorrIdx0, :)';
MovData = mesh2.vertices( CorrIdx1, :)';

C = nchoosek( 1:size(MovData, 2), 2 );
Md = createns( RefData');
tic
EvalNum = [];
Tf = zeros(3, 4, size(C, 1) );
parfor i = 1 : 1 : size(C, 1)
        [R ,T, tmpError] = RecoverRT( MovData( :, C(i, :) ), RefData( :, C(i, :) ) );
        InitData = bsxfun( @plus, R * MovData, T );
        [tmpIdx, ~] = knnsearch( Md, InitData', 'k', 1 );
        tmp = InitData - RefData(:, tmpIdx );
        tmpDist = tmp(1, :).^2 + tmp(2, :).^2 + tmp(3, :).^2;
        EvalNum = [ EvalNum length( find( tmpDist < 0.01 * 0.01 ) )];
        Tf(:, :, i) = [R T];
end
toc
[maxEvalNum Idx] = max(EvalNum);
R0 = Tf(:, 1:end-1, Idx);
T0 = Tf(:, end, Idx );
figure; 
view(3);
hold on; 
grid on; 
showPointCloud( pcData0', 'g' );
showPointCloud( pcData1', 'r' );
AftData = Loc2Glo( pcData1, R0', T0 );
showPointCloud( AftData', 'b' );
