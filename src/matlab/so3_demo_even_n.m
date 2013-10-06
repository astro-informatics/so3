% so3_demo_even_n - Run optimisation demo for even n.
%
% Simple demo to compute inverse and forward transform of complex scalar
% function, using only flmn of even n.
%
% Default usage is given by
%
%   so3_demo_even_n
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013 Martin Büttner and Jason McEwen
% See LICENSE.txt for license details

clear all;

% Define parameters
L = 64;
N = 3;

% Generate random flmn of complex signal
flmn = zeros((2*N-1)*L*L, 1);
for n = -N+1:N-1,
    if mod(n,2) == 1,
        continue
    end
    for el = abs(n):L-1,
        for m = -el:el,
            sample = rand + 1i*rand;
            sample = 2*(sample - (1+1i)/2);

            ind = so3_elmn2ind(el,m,n,L,N);
            flmn(ind) = sample;
        end
    end
end

% Compute inverse then forward transform.
f = so3_inverse(flmn, L, N, 'NMode', 'Even');
flmn_syn = so3_forward(f, L, N, 'NMode', 'Even');

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

