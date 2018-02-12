clear; clc;
addpath('Funs');
addpath('Data');
%% Set Para
opts.Win_Size = 12; % used for splitting tasks
opts.lambda1 = 1e2; % parameter before the fusion term
opts.lambda2 = 1e1; % parameter before sparse l1 norm
opts.Training_Period = [1968, 2016];
%% Get MTLData
fprintf('Getting MTL Data...\n');
File_name = strcat('Data/MTLData_win', num2str(opts.Win_Size), '_', ...
    num2str(opts.Training_Period(1)), '_', num2str(opts.Training_Period(2)), '.mat');
if ~exist( File_name, 'file')
    MTLData = GetMTLData(opts);
    save( File_name, '-v7.3', 'MTLData' );
else
    load( File_name );
end
%% Remove residues mutates less than 20% in all the tasks
Remove_pos = ReduceResidues(MTLData.TrainData, 0.2);
Remove_pos = union(Remove_pos, [2 3 6 326]);
Left_pos = setdiff( (1:329), Remove_pos );
%% Load MTLData and begin to train
TrainData = MTLData.TrainData;
Task_Num = length(TrainData);
Xmtl = cell(1, Task_Num);
Ymtl = cell(1, Task_Num);
for i = 1:Task_Num
    TrainData{i}.X = TrainData{i}.X(:, Left_pos);
    Xmtl{i} = TrainData{i}.X;
    Ymtl{i} = TrainData{i}.Y;
end
save Reduced_TrainData.mat -v7.3 TrainData Left_pos
%% Train
d = size( Xmtl{1}, 2 );
if ~exist('Weight.mat', 'file')
    %% Define graph
    [C, m_t, m_p, m_d] = DefineGraph(TrainData);
    %% Run Temporal MTL
    fprintf('Running MTL...\n');
    W0 = zeros(d, Task_Num);
    wl2 = ones(1, Task_Num);
    W_MTL = TemporalMTL_graph(Xmtl, Ymtl, C, W0, opts.lambda1, opts.lambda2, wl2);
    Site_changed = Site_Change_Statistics(Xmtl);
    save Weight.mat -v7.3 W_MTL Left_pos Site_changed
else
    [C, m_t, m_p, m_d] = DefineGraph(TrainData);
    load Weight.mat
end
%% Separate MTL and lasso
if ~exist('Weight_Compare.mat', 'file')
    W_MTL_t = TemporalMTL(Xmtl(1:m_t), Ymtl(1:m_t), zeros(d, m_t), opts.lambda1, opts.lambda2, ones(1, m_t));
    W_MTL_p = TemporalMTL(Xmtl(m_t+1:end-1), Ymtl(m_t+1:end-1), zeros(d, m_p), opts.lambda1, opts.lambda2, ones(1, m_p));
    W_MTL_d = lasso(Xmtl{end}, Ymtl{end}, 'Lambda', 0.1);
    W_lasso = zeros(d, Task_Num);
    for i = 1:Task_Num
        W_lasso(:, i) = lasso(Xmtl{i}, Ymtl{i}, 'Lambda', 0.1);
    end
    save Weight_Compare.mat -v7.3 W_MTL W_MTL_t W_MTL_p W_lasso Left_pos
else
    load Weight_Compare.mat
end
%% Display
RMSE_lasso = zeros(1, Task_Num);
Relative_Err_lasso = zeros(1, Task_Num);
for i = 1:Task_Num
    RMSE_lasso(i) = RMSE_Test(W_lasso(:, i), Xmtl{i}, Ymtl{i});
    Relative_Err_lasso(i) = RMSE_lasso(i) / mean(Ymtl{i});
end
[RMSE, Relative_Err] = Main_Test(W_MTL, TrainData);
fprintf('RMSE = %f, Relative_Err = %f\n', mean(RMSE), mean(Relative_Err));
[RMSE_t, Relative_Err_t] = Main_Test(W_MTL_t, TrainData(1:m_t));
[RMSE_p, Relative_Err_p] = Main_Test(W_MTL_p, TrainData(m_t+1:end-1));
[RMSE_d, Relative_Err_d] = Main_Test(W_MTL_d, TrainData(end));
fprintf('RMSE_separate = %f, Relative_Err = %f\n', mean([RMSE_t RMSE_p RMSE_d]), mean([Relative_Err_t Relative_Err_p Relative_Err_d]));
fprintf('RMSE_lasso = %f, Relative_Err = %f\n', mean(RMSE_lasso), mean(Relative_Err_lasso));


