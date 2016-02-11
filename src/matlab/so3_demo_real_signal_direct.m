% so3_demo_basic - Run basic demo.
%
% Simple demo to compute inverse and forward transform of complex scalar
% function, using simplest interface with default options.
%
% Default usage is given by
%
%   so3_demo_basic
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013 Martin Büttner and Jason McEwen
% See LICENSE.txt for license details

clear all;

% Define parameters
L = 16;
N = L; 

% Generate random flmn of real signal
flmn = zeros(N*L*L, 1);
% Handle n = 0
for el = 0:L-1,
% Fill fl00 with random real values
sample = rand;
sample = 2*sample - 1;

ind = so3_elmn2ind(el,0,0,L,N,'Reality',true);
flmn(ind) = sample;

% Fill fl+-m0 with conjugated random values
for m = 1:el,
sample = rand + i*rand;
sample = 2*sample - 1 - i;
ind = so3_elmn2ind(el,m,0,L,N,'Reality',true);
flmn(ind) = sample;
ind = so3_elmn2ind(el,-m,0,L,N,'Reality',true);
flmn(ind) = (-1)^m * conj(sample);
end
end

% Fill the remaining coefficients randomly
for n = 1:N-1,
for el = n:L-1,
for m = -el:el,
sample = rand + i*rand;
sample = 2*sample - 1 - i;
ind = so3_elmn2ind(el,m,n,L,N,'Reality',true);
flmn(ind) = sample;
end
end
end

% Compute inverse then forward transform.
f = so3_inverse_direct(flmn, L, N, 'Reality', true);
flmn_syn = so3_forward_direct(f, L, N, 'Reality', true);

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

