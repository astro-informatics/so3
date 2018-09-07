function [F_ssht, I_ssht] = ssht_explicit_matrix(L, method)
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
 
ssht3I = @(flm) ssht_inverse(flm, L, 'Method', method, ...
    'Reality', real);
ssht3F = @(f) ssht_inverse_adjoint(f, L, 'Method', method, ...
    'Reality', real);
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, ...
    'Grid', true, 'Method', 'MW');
 
ntheta
nphi
 
% Construct matrix F implementing forward transform.
ind = 1;
f = zeros(ntheta,nphi);
for k=1:nphi
    k
    for j=1:ntheta
        f(j,k) = 1.0;
        if (real)
          flm = ssht3F(f);
        else
          flm = ssht3F(complex(f));
        end
        F_ssht(:,ind) = flm;
        f(j,k) = 0.0;
        ind = ind+1;
    end
end
 
% Construct matrix I implementing inverse transform.
flm  = zeros(L*L, 1);
    
ind = 1;
for el = 0:L-1,
    for m = -el:el,
        ind
        ind = ssht_elm2ind(el,m)
        flm(ind) = 1.0;
        f = ssht3I(complex(flm));
        I_ssht(:,ind) = f(:);
        flm(ind) = 0.0;
        ind = ind + 1;
    end
end
 
 
% Save operators.
%filename_out = sprintf('so3_matrices_L%4.4d_N%4.4d_method%s.mat', ...
%                       L, N, method);
%save(filename_out, 'L', 'N', 'F_so3', 'I_so3');