function W = MTL_Lasso(X, Y, lambda)
Task_Num = length(X);
d = size( X{1}, 2 );
W = zeros(d, Task_Num);
for i = 1:Task_Num
    fprintf('Lasso for Task %d...\n', i);
    W(:, i) = lasso(X{i}, Y{i}, 'Lambda', lambda);
end
