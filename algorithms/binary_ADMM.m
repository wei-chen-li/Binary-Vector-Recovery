% -------------------------------------------------------------------------
% x_hat = binary_ADMM(Phi, y, lambda)
%
% Convex optimization using ADMM for binary vector recovery
% -------------------------------------------------------------------------
% Phi: M×N array
% y:   M×1 array
% x_hat: N×1 vector
%--------------------------------------------------------------------------
function x = binary_ADMM(Phi, y, lambda, varargin)

p = inputParser;
addParameter(p, 'MaxIterations', 100, @isscalar)
addParameter(p, 'x_init', [], @isvector)
parse(p, varargin{:});

max_iters = p.Results.MaxIterations;
x_init = p.Results.x_init;
clear p

[M,N] = size(Phi);

if isempty(x_init)
    x = 0.5 * ones(N,1);
else
    x = x_init;
end
z = x;
nu = zeros(N,1);

rho = 0.1;

Gram_Phi = Phi' * Phi;
I = speye(N);
A = 2 * Gram_Phi + rho * I;
A = decomposition(A, 'chol');

for iter = 1:max_iters
    x_prev = x;

    x = A \ (2 * Phi' * y + rho * z - nu);
    z = clip(x + 1/rho * nu - lambda/rho, 0,1);
    nu = nu + rho * (x - z);

    r_dual = 2 * (Gram_Phi * x - Phi' * y) + nu;
    r_prim = x - z;
    c = sqrt( (norm(r_prim,Inf) / max([norm(x,Inf) norm(z,Inf)])) / ...
              (norm(r_dual,Inf) / max([norm(2*Gram_Phi*x,Inf) norm(2*Phi'*y,Inf) norm(nu,Inf)])) );
    if c < 0.2 || c > 5
        rho = c * rho;
        A = 2 * Gram_Phi + rho * I;
        A = decomposition(A, 'chol');
    end

    if norm(x - x_prev, Inf) < 1e-6
        break
    end
end

end