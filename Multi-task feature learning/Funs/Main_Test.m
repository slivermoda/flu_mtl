function [RMSE, Relative_Err] = Main_Test(W_MTL, MTLTestData)
Task_num = length(MTLTestData);
RMSE = zeros(1, Task_num);
Relative_Err = zeros(1, Task_num);
for i = 1:Task_num
    RMSE(i) = RMSE_Test(W_MTL(:, i), MTLTestData{i}.X, MTLTestData{i}.Y);
    Relative_Err(i) = RMSE(i) / mean(MTLTestData{i}.Y);
end
