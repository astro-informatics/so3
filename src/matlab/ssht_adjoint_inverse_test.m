function [Ax, AxE] = ssht_adjoint_inverse_test(L)
% Define transforms.
real = true % this must be the case
 
sshtI = @(flm) ssht_inverse(flm, L, 'Method', 'MW', ...
    'Reality', real);
sshtIA = @(f) ssht_inverse_adjoint(f, L, 'Method', 'MW', ...
    'Reality', real);
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, ...
    'Grid', true, 'Method', 'MW');
 
ntheta
nphi

 
% Construct matrix I implementing inverse transform.
flm  = zeros(L*L, 1);
    
ind = 1;
for el = 0:L-1,
    for m = -el:el,
        ind
        ind = ssht_elm2ind(el,m)
        flm(ind) = 1.0;
        f = sshtI(complex(flm));
        I_ssht(:,ind) = f(:);
        flm(ind) = 0.0;
        ind = ind + 1;
    end
end
 

% Construct matrix F implementing forward transform.
ind = 1;
f = zeros(ntheta,nphi);
for k=1:nphi
    k
    for j=1:ntheta
        f(j,k) = 1.0;
        if (real)
          flm = sshtIA(f);
        else
          flm = sshtIA(complex(f));
        end
        F_ssht(:,ind) = flm;
        f(j,k) = 0.0;
        ind = ind+1;
    end
end

% Create random f vector
if (real)
    f = randn(ntheta,nphi);
else
    f = randn(ntheta,nphi) + 1i*randn(ntheta,nphi);
end
% Create random flmn vector
if (real)
    g = randn(ntheta,nphi);
    flm = ssht_forward(g, L, 'Method', 'MW', 'Reality', true);
else
    flm = randn(L*L, 1) + 1i*randn(L*L, 1);
end

% Compute scalar via standard transform:

Ax = sshtI(flm);

yAx = f(:)' * Ax(:)

% Compute scalar via adjoint transform:

Ay = sshtIA(f)

yAxA = Ay(:)' * flm(:)

abs(yAx - yAxA)

AyE = I_ssht' * f(:)

max(abs(Ay(:) - AyE(:)))

AxE = F_ssht' * flm(:);

max(abs(Ax(:) - AxE(:)))