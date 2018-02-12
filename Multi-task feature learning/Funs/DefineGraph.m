function [C, m_t, m_p, m_d] = DefineGraph(TrainData)
m = length(TrainData);
m_t = 0;
m_p = 0;
m_d = 0;
for i = 1:m
    if strcmp( TrainData{i}.TPYE, 'turkey' ) == 1
        m_t = m_t + 1;
    elseif strcmp( TrainData{i}.TPYE, 'pig' ) == 1
        m_p = m_p + 1;
    else
        m_d = m_d + 1;
    end
end
%% Temporal chain
C = zeros([m-1 m]);
for i = 1:m-1
    if i == m_t % Jump the gap between turkey and pig
        continue;
    end
    C(i, i) = 1;
    C(i, i+1)= -1;
end
C(m_t, :) = [];
r = size(C, 1);
%% Data similarity
for i = 1:m_t
    for j = m_t+1:m
        if TrainData{i}.A == TrainData{j}.A && TrainData{i}.B == TrainData{j}.B
            r = r + 1;
            C(r, i) = 1;
            C(r, j) = -1;
        end
    end
end
%% Drug attached at last
r = r + 1;
C(r, m - 1) = 1;
C(r, m) = -1;




