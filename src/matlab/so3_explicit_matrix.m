function [F_so3, I_so3] = so3_explicit_matrix(L, N, method)
% nsht_matrices
%
% Generate and save explicit matrices implementing forward and
% inverse so3 transforms.
%
% Default usage is given by
%
%   recsph_so3_matrices(L, N, method)
%
% where L is the harmonic band-limit, method is the shht sampling
% method.
%
% Author: Christopher Wallis (www.christophergrwallis.org)
 
% Define transforms.
real = true % this must be the case
 
so3I = @(flmn) so3_forward_adjoint_direct(flmn, L, N, 'Sampling', method, ...
    'Order', 'NegativeFirst', 'Reality', real);
so3F = @(f) so3_forward_direct(f, L, N, 'Sampling', method, ...
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
                  flmn = so3F(f);
              else
                  flmn = so3F(complex(f));
              end
              F_so3(:,ind) = flmn;
              f(i,j,k) = 0.0;
              ind = ind+1;
        end
    end
end
 
% Construct matrix I implementing inverse transform.
if (real)
    flmn  = zeros((N)*L*L, 1);
    start = 0;
else
    flmn = zeros((2*N-1)*L*L, 1);
    start = -N+1;
end

ind = 1;
for n = start:N-1,
    n
    for el = abs(n):L-1,
        for m = -el:el,
            ind
            ind = so3_elmn2ind(el,m,n,L,N, 'Order', 'NegativeFirst', 'Reality', real)
            flmn(ind) = 1.0;
            f = so3I(complex(flmn));
            I_so3(:,ind) = f(:);
            flmn(ind) = 0.0;
            ind = ind + 1;
        end
    end
end
 
 
% Save operators.
%filename_out = sprintf('so3_matrices_L%4.4d_N%4.4d_method%s.mat', ...
%                       L, N, method);
%save(filename_out, 'L', 'N', 'F_so3', 'I_so3');
