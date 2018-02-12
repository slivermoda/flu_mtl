function Best_lambda = CV_Lasso(Win_Size, Test_Year)
Lambda_Vec = 10.^(-4:2);
Method.Lasso = 1;
r = 1;
RES = cell(1);
for lambda = Lambda_Vec
    [res, ~, W] = Main_Run_OneSam(Win_Size, Test_Year, 0, 0, lambda, Method);
    RES{r} = res.Lasso;
    RES{r}.lambda = lambda;
    RES{r}.Win_Size = Win_Size;
    RES{r}.Test_Year = Test_Year;
    RES{r}.W = W;
    r = r + 1;
end
%%
minRMSE = inf;
for r = 1:length(RES)
    if minRMSE > RES{r}.RMSE
        minRMSE = RES{r}.RMSE;
        minr = r;
    end
end
Best_lambda = RES{minr}.lambda;
fprintf('minr = %d, Best_lambda = %f\n', minr, Best_lambda);



