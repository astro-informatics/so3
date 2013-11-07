% so3_demo_lower_band_limit - Run lower band-limit demo.
%
% Simple demo to compute inverse and forward transform of complex scalar
% function, using only flmn for l >= L0.
%
% Default usage is given by
%
%   so3_demo_storage_order
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013 Martin Büttner and Jason McEwen
% See LICENSE.txt for license details

clear all;

% Define parameters
L  = 64;
L0 = 32;
N  = 3;

% Generate random flmn of complex signal
flmn = zeros((2*N-1)*L*L, 1);
for n = -N+1:N-1,
    for el = max(L0, abs(n)):L-1,
        for m = -el:el,
            sample = rand + 1i*rand;
            sample = 2*(sample - (1+1i)/2);

            ind = so3_elmn2ind(el,m,n,L,N);
            flmn(ind) = sample;
        end
    end
end

% Compute inverse then forward transform.
f = so3_inverse(flmn, L, N, 'L0', L0);
flmn_syn = so3_forward(f, L, N, 'L0', L0);

% Compute maximum error in harmonic space.
maxerr = max(abs(flmn_syn - flmn))

% Compute sampling grids.
[alphas, betas, gammas, n, nalpha, nbeta, ngamma] = so3_sampling(L, N, 'Grid', true);

% Plot function on SO(3) by plotting a sphere for each value of gamma
figure(1);
for i = 1:ngamma,
    f_g = squeeze(f(i,:,:));
    subplot(1,ngamma,i)
    ssht_plot_sphere(abs(f_g), L);
    title(sprintf('g = %d; |f|', i-1))
end

