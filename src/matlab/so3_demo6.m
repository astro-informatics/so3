% so3_demo6 - Run demo6
%
% Demo to compare theoretical covariance of signal with empirical data from
% using our transform functions.
%
% Default usage is given by
%
%   so3_demo6
%
% Authors: Martin Buettner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013 Martin Buettner and Jason McEwen
% See LICENSE.txt for license details

clear all;

% Define parameters
L = 32;
N = 3;

% Generate normally distributed random flmn of complex signal
% with mean 0 and variance 1
stream = RandStream.getGlobalStream;
reset(stream);
flmn = zeros((2*N-1)*(3*L^2-N*(N-1))/3, 1);
flmn = (randn(size(flmn)) + 1i*randn(size(flmn)))/sqrt(2);

var_flmn = var(flmn);

% Compute theoretical variance of signal.
var_f_theory = 0;%zeros(1, L);
dl = zeros(2*L-1, 2*L-1);
for l = 0:L-1,
    scale = ((2*l+1)/(8*pi^2))^2;
    dl = ssht_dl(dl, L, l, pi/2);
    for m = -l:l,
        maxn = min(N-1,l);
        for n = -maxn:maxn,
            %for b = 0:L-1,
                dlmn = 0;
                for mp = -l:l, % m'
                   % We leave out the exponential contributions in alpha
                   % and gamma to Dlmn, because we are only interested in
                   % the squared modulus of Dlmn.
                   % Similarly we leave out the phase of the dlmn.
                   dlmn = dlmn + dl(L+mp,L+m)*dl(L+mp,L+n);%*exp(1i*mp*pi*(2*b+1)/(2*L-1));
                end
                %var_f_theory(b+1) = var_f_theory(b+1) + scale*abs(dlmn)^2;
                var_f_theory = var_f_theory + scale*abs(dlmn)^2;
            %end
        end
    end
end

var_f_theory = var_f_theory .* var_flmn

% Compute inverse then forward transform.
f = so3_inverse(flmn, L, N, 'Storage', 'Compact');
%f = permute(f, [2 3 1]);
%f = reshape(f, [L, (2*L-1)*(2*N-1)]).';

%var_f_data = var(f)
var_f_data = var(f(:))

max_err = max(abs(var_f_theory - var_f_data))