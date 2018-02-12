function Site_changed = Site_Change_Statistics(Xmtl)
Task_Num = length(Xmtl);
Site_changed = zeros(size(Xmtl{1}, 2), Task_Num);
for i = 1:Task_Num
    n = size(Xmtl{i}, 1);
    B = (Xmtl{i} ~= 0);
    Site_changed(:, i) = ( sum(B) / n )';
end
