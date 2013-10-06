% so3_demo_covariance - Run covariance demo.
%
% Demo to compare theoretical covariance of signal with empirical data from
% using our transform functions. The empirical covariance is computed for 
% several sets of harmonic coefficients, and the theoretical covariance is
% compared to the mean of those calculations in units of its standard
% deviation.
%
% Default usage is given by
%
%   so3_demo_covariance
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013 Martin Buettner and Jason McEwen
% See LICENSE.txt for license details

clear all;

% Define parameters
L = 16;
N = 16; % Necessary for simplification of the covariance calculation. See
        % proof further down.

var_flmn = 1; % Should we use the actual variance var(flmn) of each
              % realization here instead?

% Compute theoretical covariance of signal.
% The covariance <f(rho)f*(rho')> is 0 when rho != rho' and given by the
% following expression otherwise:
% 
% sigma^2 Sum(l,m,n) ((2*l+1)/(8pi^2))^2 |Dlmn(rho)|^2
%
% Where sigma^2 is the variance of the harmonic coefficients and Dlmn is
% the Wigner function:
%
% Dlmn(alpha, beta, gamma) = exp(-i*m*alpha)*exp(-i*n*gamma)*dlmn(beta)
%
% Note that the exponentials drop out when taking the modulus sqaured. dlmn
% are the real polar d-functions given by, which we can compute using the
% ssht_dl function. The covariance reduces to
%
%   sigma^2 Sum(l,m,n) ((2*l+1)/(8pi^2))^2 |dlmn(beta)|^2
% = sigma^2 Sum(l) (((2*l+1)/(8pi^2))^2 Sum(m,n) |dlmn(beta)|^2)
%
% The inner sum can be simplified further as follows:
%
%   Sum(m,n) |dlmn(beta)|^2
% = Sum(m,n) dlmn(beta) * dlmn(beta)*
% = Sum(m,n) dlmn(beta) * dlnm(-beta)
% = Sum(m) dlmm(beta - beta)
% = Sum(m) dlmm(0)
%
% The penultimate equality follows from the addition theorem found on pages
% 87-88 of Varshalovich. Note that this step requires the sum over n to go
% from -l to l. Hence, we can only use this simplification if N = L.
% Therefore, we obtain the following equation for the covariance:
%
% sigma^2 Sum(l,m) (((2*l+1)/(8pi^2))^2 dlmm(0)

covar_f_theory = 0; %zeros(1, L);
dl = zeros(2*L-1, 2*L-1);
for l = 0:L-1,
    scale = ((2*l+1)/(8*pi^2))^2;
    dl = ssht_dl(dl, L, l, 0);
    for m = -l:l,
        %maxn = min(N-1,l);
        %for n = -maxn:maxn,
            %for b = 0:L-1,
                %dlmn = 0;
                %for mp = -l:l, % m'
                   % We leave out the exponential contributions in alpha
                   % and gamma to Dlmn, because we are only interested in
                   % the squared modulus of Dlmn.
                   % Similarly we leave out the phase of the dlmn.
                   %dlmn = dlmn + dl(L+mp,L+m)*dl(L+mp,L+n)*exp(1i*mp*pi*(2*b+1)/(2*L-1));
                %end
                %covar_f_theory(b+1) = covar_f_theory(b+1) + scale*abs(dlmn)^2;
            %end
        %end
        covar_f_theory = covar_f_theory + scale*dl(L+m,L+m);
    end
end

covar_f_theory = covar_f_theory .* var_flmn;

runs = 100;
covar_f_data = zeros(runs,1);
    
for i = 1:runs,
    % Generate normally distributed random flmn of complex signal
    % with mean 0 and variance 1
    %stream = RandStream.getGlobalStream;
    %reset(stream);
    flmn = zeros((2*N-1)*(3*L^2-N*(N-1))/3, 1);
    flmn = (randn(size(flmn)) + 1i*randn(size(flmn)))/sqrt(2);

    % Compute inverse then forward transform.
    f = so3_inverse(flmn, L, N, 'Storage', 'Compact');
    % Reshape f, so that var(f) will give variance for each beta.
    %f = permute(f, [2 3 1]);
    %f = reshape(f, [L, (2*L-1)*(2*N-1)]).';

    %covar_f_data(i,:) = var(f)
    covar_f_data(i) = var(f(:));
end
mean_covar_f_data = mean(covar_f_data);
std_covar_f_data = std(covar_f_data);
error_in_std = abs(mean_covar_f_data - covar_f_theory)/std_covar_f_data