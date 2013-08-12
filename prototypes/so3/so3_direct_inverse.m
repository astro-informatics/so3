function [ f ] = so3_direct_inverse( L, M, N, flmn )
%SO3_DIRECT_INVERSE Inverse fast Wigner transform.
%
%   This computes the inverse Wigner transform of a function defined on the
%   rotation group SO(3). This implementation makes use of MATLAB's ifft()
%   functions.
%   
%   L    ... band limit in beta   (|b| < L)
%   M    ... band limit in alpha  (|a| < M)
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

if L < M || L < N
    error('Parameters M and N must not be greater than L.')
end
if sum(size(flmn) ~= [(2*N-1)*(3*L^2-N*(N-1))/3, 1]) > 0
    error('Parameter flmn has to be a column vector containing (2N-1)(L²-N(N-1)/3) coefficients (L²-n² for each n).')
end

% step 1
% Sum of flmn over l to obtain Fmnm'

Fmnm = zeros(2*M-1, 2*N-1, 2*L-1);

% layout of Fmnm:
% 1st index is m,  from -M+1 to M-1
% 2nd index is n,  from -N+1 to N-1
% 3rd index is m', from -L+1 to L-1

% create a lookup table that gives the flmn index of m = 0 for each l and n
noffset = zeros(2*N-1,L);
offset = 1;
for n = -N+1:N-1,
    for l = abs(n):L-1
        noffset(N+n,l+1) = offset + l;
        offset = offset + 2*l+1;
    end
end

dl = zeros(2*L-1, 2*L-1);
for l = 0:L-1,
    dl = ssht_dl(dl, L, l, pi/2);
    for m = -l:l,
        maxn = min(N-1,l);
        for n = -maxn:maxn,
            sign = 1i^(n-m);
            for mp = -l:l, % m'
                Fmnm(M+m,N+n,L+mp) = Fmnm(M+m,N+n,L+mp) +...
                    sign*(2*l+1)/(8*pi^2)*...
                    dl(L+mp,L+m)*dl(L+mp,L+n)*flmn(noffset(N+n,l+1)+m) *...
                    exp(1i*(mp)*pi/(2*L-1));
                    % this phase shift is necessary to account for the
                    % desired beta sample points (which are half a sample
                    % step size off 0)
            end
        end
    end
end

% step 2
% 3d-IFFT of Fmnm' over m', n and m to obtain f (result)
% First, we shift the 0-frequency component (m = n = m' = 0) to the 
% beginning of Fmnm' (as MATLAB expects it).
% Finally, we rearrange the dimensions of f, to obtain the desired layout.

f = ifftn(ifftshift(Fmnm));

% drop additional theta values and scale up
f = f(:, :, 1:L).*(2*L-1).*(2*M-1).*(2*N-1);

f = permute(f, [1 3 2]);

% layout of f:
% 1st index is a, from 0 to 2*M-2 (or 0 to M-1, and then -M+1 to -1)
% 2nd index is b, from 0 to L-1
% 3rd index is g, from 0 to 2*N-2 (or 0 to N-1, and then -N+1 to -1)

end