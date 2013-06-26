function [ flm ] = forward_sov_for_fft_phase( L, f )
%FORWARD_FOR_FFT Forward spherical harmonic transform using matlab's fft
%and a factoring of rotations approach.
%
%   This computes the forward spherical harmonic transform, synthesising 
%   the function using factoring of rotations. It also employs matlab's fft
%   to compute the DFT.
%   The discrepancy between the required transform and matlab's definition
%   is treated by a phase shift of function samples in phi and theta. As 
%   the function values are sampled at points in theta which are half a 
%   sample step off zero, there is also a phase shift in the m' to account 
%   for this spatial shift.
%
%   L   ... the band limit (maximum l is L-1)
%   f   ... the sampled function values

[nt,np] = size(f);
if nt ~= L || np ~= 2*L-1
    error('Parameter f has to contain L * (2*L-1) coefficients.')
end

% preparation
% phase-shift the function values in phi, to effectively cause a
% spatial shift in m such that the 0-frequency component will be in
% the middle.
phase = exp(-2*pi*1i/(2*L-1)*(0:2*L-2)*L); % shifting by -(L-1) instead of 
                                           % +L would have the same effect
phaseMat = repmat(phase, L, 1);

f = f.*phaseMat;

% step 1
% FFT of f over phi, to obtain Fm
% 1st index is t from 0 to 2L-2 (note that L to 2L-2 is the same as -L+1 to L-2)
% 2nd index is m, from -L+1 to L-1
fm = zeros(2*L-1,2*L-1); 

% fft(f,[],2) performs fft on every row of f
fm(1:L,:) = fft(f,[],2)./(2*L-1);

% step 2
% Extent the theta domain to the range [-L+1:-1] which is equivalent to
% [L:2L-1].
fm(L+1:end,:) = fm(L-1:-1:1,:);
for m = 1:2:L-1
    fm(L+1:end,L+m) = -fm(L+1:end,L+m);
    fm(L+1:end,L-m) = -fm(L+1:end,L-m);
end

% preparation
% now phase-shift the function values in theta, to effectively cause a
% spatial shift in m' such that the 0-frequency component will be in
% the middle.
phaseMat = repmat(phase.', 1, 2*L-1); 
fm = fm.*phaseMat;

% step 3
% FFT of Fm over theta, to obtain Fmm'
% then shift m = m' = 0 to the centre of the matrix to simplify step 4
% the samples are for theta values that are half a sample step off an
% alignment with 0 (which is what matlab expects). treating this
% discrepancy as a spatial shift by -pi/(2L-1), this can be corrected by a
% phase shift of the result by -pi*m'/(2L-1)
mp = repmat((-L+1:L-1).',1,2*L-1);
phase = exp(-1i*pi*mp/(2*L'-1));

% 1st index is m', from -L+1 to L-1
% 2nd index is m,  from -L+1 to L-1
fmm = fft(fm).*phase./(2*L-1);

% step 4
% Sum Fmm'' over m'', to obtain Gmm'
% 1st index is m', from -L+1 to L-1
% 2nd index is m,  from -L+1 to L-1
gmm = zeros(2*L-1, 2*L-1); 

for m=-L+1:L-1,
    for mp=-L+1:L-1,
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
            sum = sum + fmm(L+mpp,L+m)*w;
        end
        gmm(L+mp,L+m) = sum * 2*pi;
    end
end

% step 5
% Sum Gmm' over m' weighted by the Wigner functions, to obtain flm (result)
flm = zeros(L^2,1);

dl = zeros(2*L-1,2*L-1);
lmi = 1;
for l=0:L-1,
    dl = ssht_dl(dl, L, l, pi/2);
    for m=-l:l,
        sum = 0;
        for mp=-L+1:L-1,
            sum = sum + dl(L+mp,L+m)*dl(L+mp,L)*gmm(L+mp,L+m);
        end
        flm(lmi) = sum * 1i^m * sqrt((2*l+1)/(4*pi));
        lmi = lmi + 1;
    end
end