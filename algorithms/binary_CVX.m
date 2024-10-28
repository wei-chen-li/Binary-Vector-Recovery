% -------------------------------------------------------------------------
% x_hat = binary_CVX(Phi, y, 'lambda','lambda)
%
% Convex optimization for binary vector recovery
% -------------------------------------------------------------------------
% Phi: M×N×L array
% y:   M×1×L array
% x_star: N×1 vector
%--------------------------------------------------------------------------
function x_star = binary_CVX(Phi, y, varargin)

p = inputParser;
addParameter(p, 'lambda', 1, @isscalar)
parse(p, varargin{:});

lambda = p.Results.lambda;
clear p

[M,N,L] = size(Phi);

cvx_begin quiet
    variable x(N)
    minimize( (Phi * x - y)' * (Phi * x - y) + lambda * norm(x, 1) )
    subject to
        0 <= x <= 1
cvx_end

x_star = x;

end