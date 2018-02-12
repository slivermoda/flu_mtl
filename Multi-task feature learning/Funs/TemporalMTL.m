%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chain Model.
%% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = TemporalMTL(X, Y, W_initial, lambda1, lambda2, wl2)
%% lambda1: before the fusion term
%% lambda2: before the l1 norm
iter_num = 4000;
eps = 1e-5;
tau = 1;
L = 1e7;
mu = 0.01;
m = length(Y);
d = size(X{1}, 2);
% n = size(X{1}, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = zeros([m-1 m]);
for i=1:m-1
    C(i, i) = 1;
    C(i, i+1)= -1;
end
C = lambda1*C;
C = sparse(C);
A = zeros([m-1 d]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = W_initial;
gradient_old = zeros([d m]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:iter_num
    W_h = W;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_C_W = C * W_h';
    for i = 1:m-1
        A(i,:) = temp_C_W(i,:) / mu;
    end
    for i=1:m-1
        if norm(A(i,:),2) > 1
            A(i,:) = A(i,:) / norm(A(i,:),2);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute gradient
    gradient_f = zeros([d m]);
    for i = 1:m
        gradient_f(:,i) = X{i}'*( X{i}*W_h(:,i) - Y{i} );
    end
    G_f = gradient_f + A'*C;
    W_h_new = W_h - G_f/L;
    %% l1 projection
    for i = 1:m
        W_h_new(:,i) = sign(W_h_new(:,i)).*max(abs(W_h_new(:,i)) - wl2(i)*lambda2/L, 0);
    end
    %% Acceleration
    tau_new = 2 / ( iter + 3 );
    W = W_h_new + (1 - tau) / tau * tau_new * (W_h_new - W_h);
    %% Print
    Grad_Norm = norm(gradient_f, 'fro');
    if mod(iter, 100) == 0
        fprintf('Iter = %d/%d, Grad norm = %f\n', iter, iter_num, Grad_Norm);
    end
    if norm(gradient_old - gradient_f, 'fro') / Grad_Norm < eps
        fprintf('Converged.\n');
        break;
    end
end



