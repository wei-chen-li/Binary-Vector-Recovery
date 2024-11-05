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
addParameter(p, 'MaxIterations', 100, @isscalar)
addParameter(p, 'beta', NaN)
addParameter(p, 'a', 0.4, @isscalar)
addParameter(p, 'b', 0.4, @isscalar)
parse(p, varargin{:});

max_iters = p.Results.MaxIterations;
beta_hat = p.Results.beta;
a = p.Results.a;
b = p.Results.b;
clear p

[M,N,L] = size(Phi);
vecnormsqr_Phi = sum(Phi.^2, 1);

eta_0 = log(beta(a+1,b) / beta(a,b+1)) * ones(N,1);
eta_c_from = zeros(N,1,L);

if isnan(beta_hat), beta_hat = 1e2 * ones(1,1,L); end

for iter = 1:max_iters
    if exist('x_hat','var'), x_hat_old = x_hat; else, x_hat_old = inf; end

    for l = 1:L
        Phi_l = Phi(:,:,l);
        y_l = y(:,:,l);
    
        eta_c_to_l = eta_0 + sum(eta_c_from,3) - eta_c_from(:,:,l);
    
        eta_l = M_project(Phi_l, y_l, beta_hat(l), eta_c_to_l);
        x_hat = 1 ./ (1 + exp(-eta_l));
    
        eta_c_from(:,:,l) = clip(eta_l - eta_c_to_l, -1e12, 1e12);
    
        beta_hat(l) = M / ( norm(y_l - Phi_l * x_hat)^2 + vecnormsqr_Phi(:,:,l) * (x_hat .* (1-x_hat)) );
        beta_hat(l) = min(beta_hat(l), 1e10);
    end

    x_hat = 1 ./ ( 1 + exp(-eta_0 - sum(eta_c_from,3)) );

    if norm(x_hat - x_hat_old, inf) < 1e-8
        break
    end
end

end


function eta_l = M_project(Phi_l, y_l, beta_l, eta_c_to_l)

[M,N] = size(Phi_l);
x = 0.5 * ones(N,1);

for iter = 1:10
    x_old = x;
    for j = 1:N
        notj = [1:j-1 j+1:N];

        term1 = y_l - Phi_l(:,notj) * x(notj);

        z_j0 = -0.5 * beta_l * norm(term1             )^2;
        z_j1 = -0.5 * beta_l * norm(term1 - Phi_l(:,j))^2 + eta_c_to_l(j);
        x(j) = 1 / (1 + exp(z_j0 - z_j1));
    end

    if norm(x - x_old, inf) < 1e-3
        break
    end
end

eta_l = log(x ./ (1 - x));

end