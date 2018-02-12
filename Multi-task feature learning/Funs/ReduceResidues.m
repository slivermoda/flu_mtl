function pos = ReduceResidues(TrainData, Rate)
Task_num = length(TrainData);
remove_pos = cell(1, Task_num);
for i = 1:Task_num
    X = TrainData{i}.X;
    [~, remove_pos{i}] = Remove_Rare_Feature(X, Rate);
end
pos = remove_pos{1};
for i = 2:Task_num
    pos = intersect( remove_pos{i}, pos ); %% note this interscetion is remove set
end