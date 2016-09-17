function ssht_adjoint_forward_test(L)
% Define transforms.
real = true % this must be the case
 
sshtF = @(f) ssht_forward(f, L, 'Method', 'MW', ...
    'Reality', real);
sshtFA = @(flm) ssht_forward_adjoint(flm, L, 'Method', 'MW', ...
    'Reality', real);
[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, ...
    'Grid', true, 'Method', 'MW');
 
ntheta
nphi

% Create random f vector
if (real)
    f = randn(ntheta,nphi);
else
    f = randn(ntheta,nphi) + 1i*randn(ntheta,nphi);
end
% Create random flmn vector
if (real)
    flmn = zeros(L*L, 1);
    ind = 1;
    for el = 0:L-1
        flmn(ind+el) = randn();
        for m = 1:el
            val = randn() + 1i*randn();
            flmn(ind+el+m) = val;
            flmn(ind+el-m) = val * (-1)^m;
        end
        ind = ind + 2*el+1;
    end
else
    flmn = randn(L*L, 1) + 1i*randn(L*L, 1);
end

% Compute scalar via standard transform:

Ax = sshtF(f);

yAx = flmn(:)' * Ax(:)

% Compute scalar via adjoint transform:

Ay = sshtFA(flmn);

yAxA = Ay(:)' * f(:)

abs(yAx - yAxA)