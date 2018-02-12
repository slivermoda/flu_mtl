function MTLData = GetMTLData(opts)
Win_Size = opts.Win_Size;
Training_Period = opts.Training_Period;
%%
if ~exist('HI_MC_pig.mat', 'file') || ~exist('HI_MC_turkey.mat', 'file') || ~exist('HI_MC_new_drug.mat', 'file')
    error('You have to generate HI_MC first. Please go to MC folder.\n');
end
HI_pig = load('HI_MC_pig.mat');
HI_turkey = load('HI_MC_turkey.mat');
HI_drug = load('HI_MC_new_drug.mat');
HI_pig = HI_pig.HI_MC;
HI_turkey = HI_turkey.HI_MC;
HI_drug = HI_drug.HI_MC;
Train_Years = (Training_Period(1):Training_Period(2));
%% Training
TrainData_t = GetTrainData(Train_Years, Win_Size, HI_turkey, 'turkey');
TrainData_p = GetTrainData(Train_Years, Win_Size, HI_pig, 'pig');
TrainData_d = GetTrainData(Train_Years, Win_Size, HI_drug, 'drug');
TrainData = [TrainData_t TrainData_p TrainData_d];
MTLData.TrainData = TrainData;

function TestData = GetTestData(Test_Year, TrainData, HI)
A = TrainData{end}.A;
B = TrainData{end}.B;
[Sub_Mat_Train, Sub_Seq_Train] = GetSubData(HI, A, B);
[Sub_Mat_Test, Sub_Seq_Test] = GetSubData(HI, Test_Year, Test_Year);
[X_1, Y_1] = GetPairWiseData(Sub_Mat_Test, Sub_Seq_Test);
[X_2, Y_2] = GetInteractData(Sub_Mat_Train, Sub_Seq_Train, Sub_Mat_Test, Sub_Seq_Test);
TestData.X = [X_1; X_2];
TestData.Y = [Y_1; Y_2];

function MTLData = GetTrainData(Train_Year, Win_Size, HI, TPYE)
Start_Year = min( Train_Year );
End_Year = max( Train_Year );
if strcmp(TPYE, 'drug') == 1
    Task_Num = 1;
    Start_Year = 2012;
    End_Year = 2016;
    Win_Size = End_Year - Start_Year + 1;
else
    Task_Num = ( End_Year - Start_Year + 1 ) - Win_Size + 1;
end
MTLData = cell(1, Task_Num);
Empty_idx = zeros(1, Task_Num);
%%
for i = 1:Task_Num
    fprintf('Getting Task %d/%d\n', i, Task_Num);
    A = Start_Year + i - 1;
    B = A + Win_Size - 1;
    [Sub_Mat, Sub_Seq, Sub_YearList] = GetSubData(HI, A, B);
    if isempty(Sub_Mat) || size(Sub_Mat, 1) < 10 %% should at least have more than 10 viruses
        Empty_idx(i) = 1;
        continue;
    end
    MTLData{i}.SubHI = Sub_Mat;
    MTLData{i}.SubSeq = Sub_Seq;
    MTLData{i}.Virus_Num = size(Sub_Mat, 1);
    MTLData{i}.A = A;
    MTLData{i}.B = B;
    MTLData{i}.TPYE = TPYE;
    MTLData{i}.YearList = Sub_YearList;
    [ MTLData{i}.X, MTLData{i}.Y ] = GetPairWiseData(Sub_Mat, Sub_Seq);
end
MTLData(logical(Empty_idx)) = [];

function [X, Y] = GetLongDistData(HI, Win_Size, opts)
[pima, Mapping] = GetPIMA();
if opts.UseMDSdist
    Mat = HI.coor;
else
    Mat = HI.Matrix;
end
Seq = HI.Seq;
[m, n] = size(Mat);
l = length(Seq{1});
r = 0;
Max_r = m*(m-1)/2;
X = zeros(Max_r, l);
Y = zeros(Max_r, 1);
for i = 1:m
    fprintf('Long dist processing virus %d\n', i);
    for j = i+1:m
        if abs( HI.Virus_Year(j) - HI.Virus_Year(i) ) + 1 <= Win_Size
            continue;
        end
        r = r + 1;
        for k = 1:l
            if Seq{i}(k) == '-' || Seq{j}(k) == '-'
                X(r, k) = 0;
            else
                idx1 = find( Mapping == Seq{i}(k) );
                idx2 = find( Mapping == Seq{j}(k) );
                if idx1 ~= idx2
                    X(r,k) = pima(idx1,idx2);
                else
                    X(r,k) = 0;
                end
            end
        end
        if opts.UseMDSdist
            Y(r) = norm( Mat(i,:) - Mat(j,:), 2 );
        else
            Y(r) = norm( Mat(i,:) - Mat(j,:), 1 ) / n;
        end
    end
end
X = X(1:r, :);
Y = Y(1:r);

function [pima, Mapping] = GetPIMA()
%% get pima code
pima = [1,6,6,6,6,6,6,3,6,6,6,6,6,6,4,4,4,6,6,6,1;
		6,1,5,5,6,4,4,6,6,6,6,3,6,6,6,5,5,6,6,6,1;
		6,5,1,3,6,5,4,6,6,6,6,5,6,6,6,5,5,6,6,6,1;
		6,5,3,1,6,5,3,6,6,6,6,5,6,6,6,5,5,6,6,6,1;
		6,6,6,6,1,6,6,6,6,5,5,6,5,5,6,6,6,5,5,5,1;
		6,4,5,5,6,1,3,6,6,6,6,4,6,6,6,5,5,6,6,6,1;
		6,4,4,3,6,3,1,6,6,6,6,4,6,6,6,5,5,6,6,6,1;
		3,6,6,6,6,6,6,1,6,6,6,6,6,6,4,4,4,6,6,6,1;
		6,6,6,6,6,6,6,6,1,6,6,6,6,4,6,6,6,4,4,6,1;
		6,6,6,6,5,6,6,6,6,1,4,6,4,5,6,6,6,5,5,3,1;
		6,6,6,6,5,6,6,6,6,4,1,6,3,5,6,6,6,5,5,4,1;
		6,3,5,5,6,4,4,6,6,6,6,1,6,6,6,5,5,6,6,6,1;
		6,6,6,6,5,6,6,6,6,4,3,6,1,5,6,6,6,5,5,4,1;
		6,6,6,6,5,6,6,6,4,5,5,6,5,1,6,6,6,3,3,5,1;
		4,6,6,6,6,6,6,4,6,6,6,6,6,6,1,4,4,6,6,6,1;
		4,5,5,5,6,5,5,4,6,6,6,5,6,6,4,1,3,6,6,6,1;
		4,5,5,5,6,5,5,4,6,6,6,5,6,6,4,3,1,6,6,6,1;
		6,6,6,6,5,6,6,6,4,5,5,6,5,3,6,6,6,1,3,5,1;
		6,6,6,6,5,6,6,6,4,5,5,6,5,3,6,6,6,3,1,5,1;
		6,6,6,6,5,6,6,6,6,3,4,6,4,5,6,6,6,5,5,1,1;
		1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
pima = pima - 1;
Mapping = 'ARNDCQEGHILKMFPSTWYVX';

function [SubMat, SubSeq, SubYearList] = GetSubData(HI, Year_A, Year_B)
row_idx = intersect( find(HI.Virus_Year >= Year_A), find(HI.Virus_Year <= Year_B) );
SubYearList = HI.Virus_Year(row_idx);
SubMat = HI.Matrix( row_idx, : );
SubSeq = cell( 1, length(row_idx) );
for i = 1:length(row_idx)
    SubSeq{i} = HI.Seq{ row_idx(i) };
end

function [X, Y] = GetPairWiseData(Sub_Mat, Sub_Seq)
pima = GetPIMA();
[m, n] = size(Sub_Mat);
l = length(Sub_Seq{1});
X = zeros( m*(m-1)/2, l );
Y = zeros( m*(m-1)/2, 1 );
r = 1;
for i = 1:m
    for j = i+1:m
        X(r, :) = GetPairDiffX(l, Sub_Seq{i}, Sub_Seq{j}, pima);
        Y(r) = norm( Sub_Mat(i,:) - Sub_Mat(j,:), 1 ) / n;
        r = r + 1;
    end
end

function [X, Y] = GetInteractData(Mat_A, Seq_A, Mat_B, Seq_B)
[pima, Mapping] = GetPIMA();
[m1, n] = size(Mat_A);
m2 = size(Mat_B, 1);
l = length(Seq_A{1});
X = zeros( m1*m2, l );
Y = zeros( m1*m2, 1 );
r = 1;
for i = 1:m1
    for j = 1:m2
        for k = 1:l
            if Seq_A{i}(k) == '-' || Seq_B{j}(k) == '-'
                X(r, k) = 0;
            else
                idx1 = find( Mapping == Seq_A{i}(k) );
                idx2 = find( Mapping == Seq_B{j}(k) );
                if idx1 ~= idx2
                    X(r,k) = pima(idx1,idx2);
                else
                    X(r,k) = 0;
                end
            end
        end
        Y(r) = norm( Mat_A(i,:) - Mat_B(j,:), 1 ) / n;
        r = r + 1;
    end
end


