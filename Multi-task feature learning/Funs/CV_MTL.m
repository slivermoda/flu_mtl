function [Best_lambda1, Best_lambda2] = CV_MTL(Win_Size, Test_Year)
Lambda_Vec = 10.^(-4:2);
Method.MTL = 1;
r = 1;
RES = cell(1);
for lambda1 = Lambda_Vec
    for lambda2 = Lambda_Vec
        [res, W, ~] = Main_Run_OneSam(Win_Size, Test_Year, lambda1, lambda2, 0, Method);
        RES{r} = res.TempMTL;
        RES{r}.lambda1 = lambda1;
        RES{r}.lambda2 = lambda2;
        RES{r}.Win_Size = Win_Size;
        RES{r}.Test_Year = Test_Year;
        RES{r}.W = W;
        r = r + 1;
    end
end
%%
minRMSE = inf;
for r = 1:length(RES)
    if minRMSE > RES{r}.RMSE
        minRMSE = RES{r}.RMSE;
        minr = r;
    end
end
Best_lambda1 = RES{minr}.lambda1;
Best_lambda2 = RES{minr}.lambda2;
fprintf('minr = %d, Best_lambda1 = %f, Best_lambda2 = %f\n', minr, Best_lambda1, Best_lambda2);


