function so3_adjoint_forward_test(L, N)
% Define transforms.
real = true% this must be the case
nmode = 'Even' 

so3F = @(f) so3_forward_direct(f, L, N, 'Sampling', 'MW', ...
    'Order', 'NegativeFirst', 'Reality', real,  'NMode', nmode);
so3FA = @(flmn) so3_forward_adjoint_direct(flmn, L, N, 'Sampling', 'MW', ...
    'Order', 'NegativeFirst', 'Reality', real,  'NMode', nmode);
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

Ax = so3F(f);

if (real)
    upflmn = so3_unpack_flmn(flmn,L,N);
    upAx   = so3_unpack_flmn(Ax,L,N);

    yAx = upflmn(:)' * upAx(:);
else
    yAx = flmn(:)' * Ax(:);
end

%yAx = flmn(:)' * Ax(:)

% Compute scalar via adjoint transform:

Ay = so3FA(flmn);

flmn(1:500:end);
Ay(1:1000:end)

yAxA = Ay(:)' * f(:);

abs(yAx - yAxA)
