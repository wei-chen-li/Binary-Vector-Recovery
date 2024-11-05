% -------------------------------------------------------------------------
% x_hat = binary_CVX(Phi, y, 'beta',beta)
%
% Convex optimization for binary vector recovery
% -------------------------------------------------------------------------
% Phi: M×N×L array
% y:   M×1×L array
% x_star: N×1 vector
%--------------------------------------------------------------------------
function x_star = binary_CVX(Phi, y, varargin)

p = inputParser;
addParameter(p, 'MaxIterations', 100, @isscalar)
addParameter(p, 'beta', NaN)
parse(p, varargin{:});

max_iters = p.Results.MaxIterations;
beta = p.Results.beta;
clear p

[M,N,L] = size(Phi);

if L == 1 && isnan(beta)
    cvx_begin quiet
        variable x(N)
        minimize( norm(x,1) )
        subject to
            Phi * x == y
            0 <= x <= 1
    cvx_end
else
    if L == 1, max_iters = 1; end
    if isnan(beta), beta = 1e2 * ones(1,1,L); end

    for iter = 1:max_iters
        Phi_tilde = reshape(permute(Phi .* sqrt(beta), [1 3 2]), [M*L N]);
        y_tilde =   reshape(permute(Phi .* sqrt(beta), [1 3 2]), [M*L 1]);

        if exist('x','var'), x_old = x; else, x_old = Inf; end
        
        cvx_begin quiet
            variable x(N)
            minimize( norm(x,1) )
            subject to
                (Phi_tilde * x - y_tilde)' * (Phi_tilde * x - y_tilde) <= 2.5 * M*L;
                0 <= x <= 1
        cvx_end

        for l = 1:L
            beta(l) = M / norm(y - Phi(:,:,l) * x)^2;
            beta(l) = min(beta(l), 1e10);
        end

        if norm(x - x_old, inf) < 1e-8
            break
        end
    end
end

x_star = x;

end