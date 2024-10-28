% -------------------------------------------------------------------------
% x_hat = sprase_CVX(Phi, y, 'lambda','lambda)
%
% Convex optimization for binary vector recovery
% -------------------------------------------------------------------------
% Phi: M×N×L array
% y:   M×1×L array
% x_star: N×1 vector
%--------------------------------------------------------------------------
function x_star = sparse_CVX(Phi, y, varargin)

p = inputParser;
addParameter(p, 'lambda', 1, @isscalar)
parse(p, varargin{:});

lambda = p.Results.lambda;
clear p

[M,N] = size(Phi);

cvx_begin quiet
    variable x(N)
    minimize( (Phi * x - y)' * (Phi * x - y) + lambda * norm(x, 1) )
cvx_end

x_star = x;

end