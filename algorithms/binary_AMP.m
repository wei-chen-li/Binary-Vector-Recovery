% -------------------------------------------------------------------------
% x_hat = binary_AMP(Phi, y, 'beta',beta, 'a',a, 'b',b)
%
% Approximate message passing for binary vector recovery
% -------------------------------------------------------------------------
% Phi: M×N×L array
% y:   M×1×L array
% x_hat: N×1 vector
%--------------------------------------------------------------------------
function x_hat = binary_AMP(Phi,y, varargin)

p = inputParser;
addParameter(p, 'MaxIterations', 30, @isscalar)
addParameter(p, 'MaxIterationsInnerLoop', 100, @isscalar)
addParameter(p, 'beta', NaN)
addParameter(p, 'a', 0.4, @isscalar)
addParameter(p, 'b', 0.4, @isscalar)
parse(p, varargin{:});

max_iters = p.Results.MaxIterations;
max_iters_inner = p.Results.MaxIterationsInnerLoop;
beta_hat = p.Results.beta;
a = p.Results.a;
b = p.Results.b;
clear p

[M,N,L] = size(Phi);
vecnormsqr_Phi = sum(Phi.^2, 1);

eta_c_to_r = log(beta(a+1,b) / beta(a,b+1)) * ones(N,1);

if isnan(beta_hat), beta_hat = 1e2 * ones(1,1,L); end
beta_hat = min(beta_hat, 1e10);

for iter = 1:max_iters
    if exist('x_hat','var'), x_hat_old = x_hat; else, x_hat_old = inf; end

    x_hat = M_project(Phi, y, beta_hat, eta_c_to_r, 0.5*ones(N,1), max_iters_inner);

    for l = 1:L    
        beta_hat(l) = M / ( norm(y(:,:,l) - Phi(:,:,l) * x_hat)^2 + vecnormsqr_Phi(:,:,l) * (x_hat .* (1-x_hat)) );
        beta_hat(l) = min(beta_hat(l), 1e10);
    end

    if norm(x_hat - x_hat_old, inf) < 1e-8
        break
    end
end

end


function x_hat = M_project(Phi, y, beta_hat, eta_c_to_r, x_hat_init, max_iters_inner)

[M,N,L] = size(Phi);
beta_hat = reshape(beta_hat, 1,1,L);

x_hat = x_hat_init;

for iter = 1:max_iters_inner
    x_hat_old = x_hat;

    for j = 1:N
        notj = [1:j-1 j+1:N];

        term1 = y - pagemtimes(Phi(:,notj,:), x_hat(notj));

        z_j0 = sum(-0.5 * beta_hat .* vecnorm(term1             ).^2, 3);
        z_j1 = sum(-0.5 * beta_hat .* vecnorm(term1 - Phi(:,j,:)).^2, 3) + eta_c_to_r(1);
        x_hat(j) = 1 / (1 + exp(z_j0 - z_j1));
    end

    if norm(x_hat - x_hat_old, inf) < 1e-6
        break
    end
end

end