% so3_demo_demo_plot_wigner - Run demo to create a plot of a Wigner
% function.
%
% Plot Wigner functions on SO(3) by plotting a sphere for each sample
% of gamma.
%
% Default usage is given by
%
%   so3_demo_demo_plot_wigner
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013 Martin Büttner and Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Define parameters.
L = 64;
N = 3;
el = 4;
m  = 2;
n  = 1;
type = 'colour';

% Generate Wigner coefficients
flmn = zeros((2*N-1)*L*L, 1);
ind = so3_elmn2ind(el, m, n, L, N);
flmn(ind) = 1.0;

% Compute function on SO(3)
f = so3_inverse(complex(real(flmn), imag(flmn)), L, N);

% Compute sampling grids
[alphas, betas, gammas, n, nalpha, nbeta, ngamma] = so3_sampling(L, N, 'Grid', true);

% Plot function on SO(3) by plotting a sphere for each value of gamma
figure(1);
for i = 1:ngamma,
    f_g = squeeze(f(i,:,:));
    subplot(2,ngamma,i)
    ssht_plot_sphere(real(f_g), L, 'Type', type,...
        'ColourBar', false, 'Lighting', true);
    title(sprintf('g = %d; Re(f)', i-1))
    subplot(2,ngamma,ngamma+i)
    ssht_plot_sphere(imag(f_g), L, 'Type', type,...
        'ColourBar', false, 'Lighting', true);
    title(sprintf('g = %d; Im(f)', i-1))
end
