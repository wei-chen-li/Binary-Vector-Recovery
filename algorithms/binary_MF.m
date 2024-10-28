% -------------------------------------------------------------------------
% x_hat = binary_MF(Phi, y, 'a',a, 'b',b)
%
% Mean field approximate inference for binary vector recovery
% -------------------------------------------------------------------------
% Phi: M×N×L array
% y:   M×1×L array
% x_hat: N×1 vector
%--------------------------------------------------------------------------
function x_hat = binary_MF(Phi, y, varargin)

p = inputParser;
addParameter(p, 'MaxIterations', 100, @isscalar)
addParameter(p, 'a', 0.4, @isscalar)
addParameter(p, 'b', 0.4, @isscalar)
parse(p, varargin{:});

max_iters = p.Results.MaxIterations;
a = p.Results.a;
b = p.Results.b;
clear p

[M,N,L] = size(Phi);
vecnormsqr_Phi = sum(Phi.^2, 1);

beta_hat = zeros(1,1,L);
x_hat = 0.5 * ones(N,1);

for iter = 1:max_iters
    for l = 1:L
        y_l = y(:,:,l);
        Phi_l = Phi(:,:,l);
        beta_hat(l) = M / ( norm(y_l - Phi_l * x_hat)^2 + vecnormsqr_Phi(:,:,l) * (x_hat .* (1-x_hat)) );
    end

    log_pi     = psi(    x_hat + a) - psi(a+b+1);
    log_neg_pi = psi(1 - x_hat + b) - psi(a+b+1);

    x_hat_old = x_hat;

    y_notj = y - pagemtimes(Phi, x_hat) + Phi .* x_hat';
    value = -0.5 * beta_hat .* (vecnormsqr_Phi - 2 * sum(Phi .* y_notj, 1));
    value = reshape(sum(value,3), N,1);
    x_hat = 1 ./ ( 1 + exp(log_neg_pi - value - log_pi) );

    if norm(x_hat - x_hat_old, inf) < 1e-8
        break
    end
end

end