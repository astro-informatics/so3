function [ f ] = so3_via_ssht_inverse( L, N, flmn )
%SO3_VIA_SSHT_INVERSE Inverse fast Wigner transform exploiting the fast
%spin spherical harmonic transform.
%
%   This computes the inverse Wigner transform of a function defined on the
%   rotation group SO(3). This implementation makes use of MATLAB's ifft()
%   functions.
%   
%   L    ... band limit in beta   (|b| < L)
%            this implicitly imposes a band limit M = L on alpha
%   N    ... band limit in gamma  (|c| < N)
%   flmn ... the coefficients of the Wigner D-function conjugates.
%            To save memory and allow for easy exploitation of the fast
%            SSHT, the coefficients are expected in an unrolled form as
%            follows:
%            The flmn array consists of 2*N-1 sections, each corresponding
%            to a value of n from -N+1 to N+1.
%            Each of those sections contains L-|n| subsections, each
%            corresponding to a value of l from |n| to L-1.
%            Each of those subsections lists 2*l+1 values of m ranging 
%            from -l+1 to l-1.
%            Note that each n-section contains L²-n² values.
%            This sums to a total of (2N-1)(L²-N(N-1)/3) values.

if L < N
    error('Band limit N must not be greater than L.')
end

if sum(size(flmn) ~= [(2*N-1)*(3*L^2-N*(N-1))/3, 1]) > 0
    error('Parameter flmn has to be a column vector containing (2N-1)(L²-N(N-1)/3) coefficients (L²-n² for each n).')
end

% step 1
% SSHT of flmn over l and m to obtain fn(a,b)

% Create a mask of scaling values that can be applied to an unrolled flm
% array.
scale = zeros(L^2,1);
offset = 1; % marks the beginning of the current l-subsection
for l = 0:L-1
    scale(offset:offset+2*l) = sqrt((2*l+1)/(16*pi^3));
    offset = offset + 2*l+1;
end


fn = zeros(2*N-1, L, 2*L-1);

offset = 1; % marks the beginning of the current n-section
for n = -N+1:N-1
    % Slice the current section from flmn and pad with zeros
    flm = [zeros(n^2,1); flmn(offset:offset+L^2-n^2-1)];
    % Scale flm
    flm = flm.*scale;
    
    fn(N+n,:,:) = (-1)^n*ssht_inverse(flm, L, 'Spin', -n);
    offset = offset + L^2 - n^2;
end

% layout of fn:
% 1st index is n, from -N+1 to N-1
% 2nd index is b, from 0 to L-1
% 3rd index is a, from 0 to 2*L-2 (or 0 to L-1, and then -L+1 to -1)

% step 2
% IFFT of fn(a,b) over n to obtain f(a,b,g) (result)
% First, we shift the 0-frequency component (n = 0) to the 
% beginning of fn (as MATLAB expects it).
% Finally, we rearrange the dimensions of f, to obtain the desired layout.
f = ifft(ifftshift(fn,1))*(2*N-1);

f = permute(f,[3 2 1]);

% layout of f:
% 1st index is a, from 0 to 2*L-2 (or 0 to L-1, and then -L+1 to -1)
% 2nd index is b, from 0 to L-1
% 3rd index is g, from 0 to 2*N-2 (or 0 to N-1, and then -N+1 to -1)

end