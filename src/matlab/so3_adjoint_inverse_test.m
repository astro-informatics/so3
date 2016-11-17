function [Ax, AxE] = so3_adjoint_inverse_test(L, N)
% Define transforms.
real = true % this must be the case
 
so3I = @(flmn) so3_inverse_direct(flmn, L, N, 'Sampling', 'MW', ...
    'Order', 'NegativeFirst', 'Reality', real);
so3IA = @(f) so3_inverse_adjoint_direct(f, L, N, 'Sampling', 'MW', ...
    'Order', 'NegativeFirst', 'Reality', real);
[alphas, betas, gammas, n, nalpha, nbeta, ngamma] = so3_sampling(L, N, ...
    'Grid', true, 'Method', 'MW');
 
nalpha
nbeta
ngamma

 
% Construct matrix F implementing forward transform.
ind = 1;
f = zeros(ngamma,nbeta,nalpha);
for k=1:nalpha
    k
    for j=1:nbeta
        for i=1:ngamma
              f(i,j,k) = 1.0;
              if (real)
                  flmn = so3IA(f);
              else
                  flmn = so3IA(complex(f));
              end
              F_so3(:,ind) = flmn;
              f(i,j,k) = 0.0;
              ind = ind+1;
        end
    end
end

% Create random f vector
if (real)
    f = randn(ngamma,nbeta,nalpha);
else
    f = randn(ngamma,nbeta,nalpha) + 1i*randn(ngamma,nbeta,nalpha);
end
% Create random flmn vector
if (real)
    g = randn(ngamma,nbeta,nalpha);
    flmn = so3_forward(g, L, N, 'Sampling', 'MW', 'Order', 'NegativeFirst', 'Reality', true);
else
    flmn = randn((2*N-1)*L*L, 1) + 1i*randn((2*N-1)*L*L, 1);
end

% Compute scalar via standard transform:

Ax = so3I(flmn);

yAx = f(:)' * Ax(:)

% Compute scalar via adjoint transform:

Ay = so3IA(f);

upflmn = so3_unpack_flmn(flmn,L,N);
upAy   = so3_unpack_flmn(Ay,L,N);

yAxA = upAy(:)' * upflmn(:)

test_dot = flmn(:)' * flmn(:)

yAx - yAxA

%AxE = F_so3' * flmn(:)

%max(abs(Ax(:) - AxE(:)))

