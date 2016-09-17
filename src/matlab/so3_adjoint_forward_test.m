function so3_adjoint_forward_test(L, N)
% Define transforms.
real = false % this must be the case
 
so3F = @(f) so3_forward_direct(f, L, N, 'Sampling', 'MW', ...
    'Order', 'NegativeFirst', 'Reality', real);
so3FA = @(flmn) so3_forward_adjoint_direct(flmn, L, N, 'Sampling', 'MW', ...
    'Order', 'NegativeFirst', 'Reality', real);
[alphas, betas, gammas, n, nalpha, nbeta, ngamma] = so3_sampling(L, N, ...
    'Grid', true, 'Method', 'MW');
 
nalpha
nbeta
ngamma

% Create random f vector
if (real)
    f = randn(ngamma,nbeta,nalpha);
else
    f = randn(ngamma,nbeta,nalpha) + 1i*randn(ngamma,nbeta,nalpha);
end
% Create random flmn vector
if (real)
    flmn = randn((N)*L*L, 1) + 1i*randn((N)*L*L, 1);
else
    flmn = randn((2*N-1)*L*L, 1) + 1i*randn((2*N-1)*L*L, 1);
end

% Compute scalar via standard transform:

Ax = so3F(f);

yAx = flmn(:)' * Ax(:)

% Compute scalar via adjoint transform:

Ay = so3FA(flmn);

yAxA = Ay(:)' * f(:)

abs(yAx - yAxA)