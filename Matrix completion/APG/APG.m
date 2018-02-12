function Y = APG(X, M, known, eps, lambda, L)
%% min  g(x)+h(x)    g(x) = 1/2 * L .* X .*X + 1/2*||X_\Omega - M_\Omega||_F^2, h(x) = lambda*||X||_*
Max_iter = 500;
Lips = 1;
lastX = X;
Y = lastX;
tlast = 1;
for k = 1:Max_iter
    Count = 1;
    while true
        [u, sigma, v] = svd(Y - 1/Lips * ( L.*Y + (Y-M).*known));
        X = u * max( sigma - lambda / Lips, 0 ) * v';
        if obj_val(X, lambda, L, M, known) <= obj_val(lastX, lambda, L, M, known)
            break;
        end
        Lips = Lips*1.1;
        Count = Count + 1;
        if Count > 100
            return;
        end
    end
    t = (1 + sqrt( 1 + 4*tlast^2 )) / 2;
    Y = X + ( tlast - 1 ) / t * (X - lastX);
    tlast = t;       
    history.objval(k) = obj_val(Y, lambda, L, M, known);
    Err = norm(Y - lastX, 'fro') / norm(lastX, 'fro');
    if mod(k, 100) == 0
%         fprintf('iter = %d, obj = %f, Err = %f\n', k, history.objval(k), Err );
    end
    if( k >= 2 && Err < eps )
        fprintf('Converged.\n');
        break;
    end
    lastX = X;
end

function obj = obj_val(X, lambda, L, M, known)
obj = lambda*trace_norm(X) + 1/2*sum(sum(L.*X.*X)) + 1/2*norm((X-M).*known,'fro')^2;

function obj = trace_norm(X)
sigma = svd(X);
obj = sum(abs(sigma));




