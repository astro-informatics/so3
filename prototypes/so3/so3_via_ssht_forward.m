function [ flmn ] = so3_via_ssht_forward( L, N, f )
%SO3_VIA_SSHT_FORWARD Forward fast Wigner transform exploiting the fast
%spin spherical harmonic transform.
%
%   This computes the forward Wigner transform of a function defined on the
%   rotation group SO(3). This implementation makes use of MATLAB's fft()
%   functions.
%   
%   L ... band limit in beta   (|b| < L)
%         this implicitly imposes a band limit M = L on alpha
%   N ... band limit in gamma  (|c| < N)
%   f ... sampled function values. The dimensions of f correspond to alpha,
%         beta, gamma and the values should correspond to sampling points
%         as given by so3_sampling().
%
%   To save memory and allow for easy exploitation of the fast SSHT, the 
%   returned coefficients flmn are expected in an unrolled form as follows:
%   The flmn array consists of 2*N-1 sections, each corresponding
%   to a value of n from -N+1 to N+1.
%   Each of those sections contains L-|n| subsections, each
%   corresponding to a value of l from |n| to L-1.
%   Each of those subsections lists 2*l+1 values of m ranging 
%   from -l+1 to l-1.
%   Note that each n-section contains L²-n² values.
%   This sums to a total of (2N-1)(L²-N(N-1)/3) values.

if L < N
    error('Parameter N must not be greater than L.')
end
[na, nb, ng] = size(f);
if na ~= 2*L-1 || nb ~= L || ng ~= 2*N-1
    error('Parameter f has to contain (2*L-1) * L * (2*N-1) coefficients.')
end

% layout of f:
% 1st index is a, from 0 to 2*M-2 (or 0 to M-1, and then -M+1 to -1)
% 2nd index is b, from 0 to L-1
% 3rd index is g, from 0 to 2*N-2 (or 0 to N-1, and then -N+1 to -1)

% step 1
% FFT of f(a,b,g) over g to obtain fn(a,b)

fn = fftshift(fft(f,[],3),3)*2*pi/(2*N-1);

% layout of fn:
% 1st index is a, from 0 to 2*L-2 (or 0 to L-1, and then -L+1 to -1)
% 2nd index is b, from 0 to L-1
% 3rd index is n, from -N+1 to N-1

% step 2
% SSHT of fn(a,b) over a and b to obtain flmn (result)

flmn = zeros((2*N-1)*(3*L^2-N*(N-1))/3,1);

% Create a mask of scaling values that can be applied to an unrolled flm
% array.
scale = zeros(L^2,1);
offset = 1; % marks the beginning of the current l-subsection
for l = 0:L-1
    scale(offset:offset+2*l) = sqrt(4*pi/(2*l+1));
    offset = offset + 2*l+1;
end

offset = 1; % marks the beginning of the current n-section
for n = -N+1:N-1
    flm = (-1)^n*ssht_forward(fn(:,:,N+n).', L, 'Spin', -n);
    % Scale flm
    flm = flm.*scale;
    
    flmn(offset:offset+L^2-n^2-1) = flm(n^2+1:end);
    
    offset = offset + L^2 - n^2;
end

end