function MakeFile(Test_Year, RES, W_MTL, W_Lasso)
filename = strcat('Res/res', num2str(Test_Year), '.txt');
fid = fopen(filename, 'w');
%%
RES.Lasso
fprintf(fid, 'MTL:\n' );
fprintf(fid, 'RMSE = %f\n', RES.TempMTL.RMSE );
fprintf(fid, 'Acc = %f\n', RES.TempMTL.Acc );
fprintf(fid, 'Sen = %f\n', RES.TempMTL.Sen );
fprintf(fid, 'Spe = %f\n', RES.TempMTL.Spe );
fprintf(fid, 'nnz = %d\n', nnz(W_MTL(:,end)) );
%%
fprintf(fid, 'Lasso:\n' );
fprintf(fid, 'RMSE = %f\n', RES.Lasso.RMSE );
fprintf(fid, 'Acc = %f\n', RES.Lasso.Acc );
fprintf(fid, 'Sen = %f\n', RES.Lasso.Sen );
fprintf(fid, 'Spe = %f\n', RES.Lasso.Spe );
fprintf(fid, 'nnz = %d\n', nnz(W_Lasso(:,end)) );
%%
filename = strcat('Res/RES', num2str(Test_Year), '.mat');
save(filename, 'RES', 'W_MTL', 'W_Lasso');