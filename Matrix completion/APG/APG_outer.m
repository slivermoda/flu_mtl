function X = APG_outer(M, known, lambda)
%%
[m, n] = size(M);
eps = 1e-3;
%% Kernal weight
L = zeros(m, n);
% H = 1e-1 / 2;
% for i = 1:m
%     for j = 1:n
%         L(i,j) = 1 / exp( - abs( i/m - j/n ) / H );
%     end
% end
% L = L*1e-10;
%%
X = M;
X = APG(X, M, known, eps, lambda, L);

