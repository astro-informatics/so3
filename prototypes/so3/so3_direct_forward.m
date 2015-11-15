function [ flmn ] = so3_direct_forward( L, M, N, f )
%SO3_DIRECT_FORWARD Forward fast Wigner transform.
%
%   This computes the forward Wigner transform of a function defined on the
%   rotation group SO(3). This implementation makes use of MATLAB's fft()
%   functions.
%   
%   L ... band limit in beta   (|b| < L)
%   M ... band limit in alpha  (|a| < M)
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

if L < M || L < N
    error('Parameters M and N must not be greater than L.')
end
[na, nb, ng] = size(f);
if na ~= 2*M-1 || nb ~= L || ng ~= 2*N-1
    error('Parameter f has to contain (2*M-1) * L * (2*N-1) coefficients.')
end

% layout of f:
% 1st index is a, from 0 to 2*M-2 (or 0 to M-1, and then -M+1 to -1)
% 2nd index is b, from 0 to L-1
% 3rd index is g, from 0 to 2*N-2 (or 0 to N-1, and then -N+1 to -1)

% step 1
% 2d-FFT of f over a and g to obtain Fmn(b)
% first, rearrange dimensions of f for use with fft2().
% new layout of f:
% 1st index is a, from 0 to 2*M-2 (or 0 to M-1, and then -M+1 to -1)
% 2nd index is g, from 0 to 2*N-2 (or 0 to N-1, and then -N+1 to -1)
% 3rd index is b, from 0 to L-1

f = permute(f, [1 3 2]);

Fmn = zeros(2*M-1, 2*N-1, 2*L-1);
Fmn(:,:,1:L) = fft2(f)/(2*M-1)/(2*N-1);

% layout of Fmn:
% 1st index is m, from 0 to 2*M-2 (or 0 to M-1, and then -M+1 to -1)
% 2nd index is n, from 0 to 2*N-2 (or 0 to N-1, and then -N+1 to -1)
% 3rd index is b, from 0 to 2*L-2 (or 0 to L-1, and then -L+1 to -1)


% step 2
% Extend the beta domain to the range [-L+1:-1] which is equivalent to
% [L:2L-1].
Fmn(:,:,L+1:end) = Fmn(:,:,L-1:-1:1);
m = (1:2:M-1);
Fmn(1+m,:,L+1:end) = -Fmn(1+m,:,L+1:end);
Fmn(2*M-m,:,L+1:end) = -Fmn(2*M-m,:,L+1:end);
n = (1:2:N-1);
Fmn(:,1+n,L+1:end) = -Fmn(:,1+n,L+1:end);
Fmn(:,2*N-n,L+1:end) = -Fmn(:,2*N-n,L+1:end);
    
% step 3
% FFT of Fmn(b) over b to obtain Fmnm'
% then shift m = n = m' = 0 to the center of the matrix to simplify steps 4
% and 5
% the samples are for beta values that are half a sample step off an
% alignment with 0 (which is what MATLAB expects). treating this
% discrepancy as a spatial shift by -pi/(2L-1), this can be corrected by a
% phase shift of the result by -pi*m'/(2L-1)
mp = repmat(shiftdim((-L+1:L-1),-1),[2*M-1, 2*N-1, 1]);
phase = exp(-1i*pi*mp/(2*L-1));
Fmnm = fftshift(fft(Fmn, [], 3)).*phase/(2*L-1);

% layout of Fmnm:
% 1st index is m,  from -M+1 to M-1
% 2nd index is n,  from -N+1 to N-1
% 3rd index is m', from -L+1 to L-1

% step 4
% Convolution of Fmnm' with w(m'' - m') to obtain Gmnm'
Gmnm = zeros(size(Fmnm));
for m = -M+1:M-1
    for n = -N+1:N-1
        for mp = -L+1:L-1
            sum = 0;
            for mpp=-L+1:L-1,
                dm = mpp-mp;
                if abs(dm) == 1 % this could be avoided by splitting up the loop
                    w = 1i*dm*pi/2;
                elseif mod(dm,2) == 0
                    w = 2/(1-dm^2);
                else
                    w = 0;
                end
                sum = sum + Fmnm(M+m, N+n, L+mpp)*w;
            end
            Gmnm(M+m, N+n, L+mp) = sum * (2*pi)^2;
        end
    end
end

% step 5
% Sum of Gmnm' over m' to obtain flmn
flmn = zeros((2*N-1)*(3*L^2-N*(N-1))/3,1);

% create a lookup table that gives the index of m = 0 for each l and n
noffset = zeros(2*N-1,L);
offset = 1;
for n = -N+1:N-1,
    for l = abs(n):L-1
        noffset(N+n,l+1) = offset + l;
        offset = offset + 2*l+1;
    end    
end

dl = zeros(2*L-1,2*L-1);
for l = 0:L-1,
    dl = ssht_dl(dl, L, l, pi/2);
    maxn = min(N-1,l);
    for n = -maxn:maxn,
        for m = -l:l,
            sum = 0;
            for mp=-L+1:L-1,
                sum = sum + dl(L+mp,L+m)*dl(L+mp,L+n)*Gmnm(M+m, N+n, L+mp);
            end
            flmn(noffset(N+n,l+1)+m) = sum * 1i^(m-n);
        end
    end
end

end