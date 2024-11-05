% -------------------------------------------------------------------------
% x_hat = sparse_Bayesian(Phi, y)
%
% Convex optimization for binary vector recovery
% -------------------------------------------------------------------------
% Phi: M×N array
% y:   M×1 array
% x_hat: N×1 vector
%--------------------------------------------------------------------------
function x_hat = sparse_Bayesian(Phi, y, varargin)

p = inputParser;
addParameter(p, 'MaxIterations',100, @isscalar)
parse(p, varargin{:});

max_iters = p.Results.MaxIterations;
clear p

[M,N] = size(Phi);

beta = 1e5;
alpha = ones(N,1);

for iter = 1:max_iters
    Sigma0 = spdiags(1./alpha, 0, N,N);

    Sigmay = 1/beta * eye(M) + Phi * Sigma0 * Phi';
    PhiSigma0 = Phi * Sigma0;
    Sigmax = Sigma0 - PhiSigma0' * (Sigmay \ PhiSigma0);
    mux = Sigmax * (Phi' * y) * beta;

    for i = 1:N
        alpha(i) = 1 / (Sigmax(i,i) + mux(i)^2);
    end

    beta = M / (norm(y - Phi * mux)^2 + trace(Phi * Sigmax * Phi'));
end

x_hat = mux;

end