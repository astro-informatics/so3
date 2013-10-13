% so3_demo3 - Run demo3
%
% Simple demo to compute inverse and forward transform of real scalar
% function, using simplest interface with default options.
%
% Default usage is given by
%
%   so3_demo3
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

% Generate random flmn of real signal
flmn = zeros((2*N-1)*L*L, 1);
for n = -N+1:N-1,
    for el = abs(n):L-1,
        for m = 0:el,
            sample = rand + 1i*rand;
            sample = 2*(sample - (1+1i)/2);

            if m == 0 && n == 0
                ind = so3_elmn2ind(el,m,n,L,N);
                flmn(ind) = real(sample);
            else
                ind = so3_elmn2ind(el,m,n,L,N);
                flmn(ind) = sample;
                ind = so3_elmn2ind(el,-m,-n,L,N);
                flmn(ind) = (-1)^(m+n) * conj(sample);
            end
        end
    end
end

% Compute inverse then forward transform.
f = so3_inverse(flmn, L, N);
flmn_syn = so3_forward(f, L, N);

% Compute maximum error in harmonic space.
maxerr = max(abs(flmn_syn - flmn))
% Compute maximum imaginary part in SO(3).
maximag = max(imag(f(:)))

% Compute sampling grids.
[alphas, betas, gammas, n, nalpha, nbeta, ngamma] = so3_sampling(L, N, 'Grid', true);

% Plot function on SO(3) by plotting a sphere for each value of gamma
figure(1);
for i = 1:ngamma,
    f_g = squeeze(f(i,:,:));
    subplot(1,ngamma,i)
    ssht_plot_sphere(real(f_g), L);
    title(sprintf('g = %d; Re(f)', i-1))
end

