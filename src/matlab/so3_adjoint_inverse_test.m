function [Ax, AxE] = so3_adjoint_inverse_test(L, N)
% Define transforms.
real = false % this must be the case
nmode = 'Even' 

so3I = @(flmn) so3_inverse_direct(flmn, L, N, 'Sampling', 'MW', ...
    'Order', 'NegativeFirst', 'Reality', real, 'Nmode', nmode);
so3IA = @(f) so3_inverse_adjoint_direct(f, L, N, 'Sampling', 'MW', ...
    'Order', 'NegativeFirst', 'Reality', real, 'Nmode', nmode);
[alphas, betas, gammas, n, nalpha, nbeta, ngamma] = so3_sampling(L, N, ...
    'Grid', true, 'Method', 'MW');
 
nalpha
nbeta
ngamma


% Create random f vector
if (real)
    f = randn(ngamma,nbeta,nalpha);
    flmn = so3_forward(f, L, N, 'Sampling', 'MW', 'Order', 'NegativeFirst', 'Reality', true, ...
         'NMode', nmode);
    f = so3_inverse(flmn, L, N, 'Sampling', 'MW', 'Order', 'NegativeFirst', 'Reality', true, ...
         'NMode', nmode);
else
    f = randn(ngamma,nbeta,nalpha) + 1i*randn(ngamma,nbeta,nalpha);
    flmn = so3_forward(f, L, N, 'Sampling', 'MW', 'Order', 'NegativeFirst', 'Reality', false, ...
         'NMode', nmode);
    f = so3_inverse(flmn, L, N, 'Sampling', 'MW', 'Order', 'NegativeFirst', 'Reality', false, ...
         'NMode', nmode);
end
% Create random flmn vector
if (real)
    g = randn(ngamma,nbeta,nalpha);
    flmn = so3_forward(g, L, N, 'Sampling', 'MW', 'Order', 'NegativeFirst', 'Reality', true);
else
    g = randn(ngamma,nbeta,nalpha) + 1i*randn(ngamma,nbeta,nalpha);
    flmn = so3_forward(g, L, N, 'Sampling', 'MW', 'Order', 'NegativeFirst', 'Reality', false, ...
         'NMode', nmode);
end


% Compute scalar via standard transform:

Ax = so3I(flmn);

yAx = f(:)' * Ax(:);

% Compute scalar via adjoint transform:

Ay = so3IA(f);

if (real)
    upflmn = so3_unpack_flmn(flmn,L,N);
    upAy   = so3_unpack_flmn(Ay,L,N);

    yAxA = upAy(:)' * upflmn(:);
else
    yAxA = Ay(:)' * flmn(:);
end

test_dot = flmn(:)' * flmn(:);

yAx - yAxA

%AxE = F_so3' * flmn(:);

%max(abs(yAx(:) - AxE(:)))
