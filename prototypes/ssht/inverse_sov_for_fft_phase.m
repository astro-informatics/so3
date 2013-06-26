function [ f ] = inverse_sov_for_fft_phase( L, flm )
%INVERSE_FOR_FFT Inverse spherical harmonic transform using matlab's ifft
%and a factoring of rotations approach.
%.
%   This computes the inverse spherical harmonic transform, synthesising 
%   the function using factoring of rotations (as outlined in equations
%   (12) to (14) of the ssht paper). It also employs matlab's ifft to 
%   compute the DFT.
%   The discrepancy between the required transform and matlab's definition
%   is treated by a phase shift of function samples in phi and theta. As 
%   the function values are sampled at points in theta which are half a 
%   sample step off zero, there is also a phase shift in the m' to account 
%   for this spatial shift.
%
%   L   ... the band limit (maximum l is L-1)
%   flm ... the coefficients of the spherical harmonics
%           to avoid wasting memory, these should be supplied in an
%           unrolled array of the following format:
%           [(0,0) (1,-1) (1,0) (1,1) (2,-2) (2,-1) ... ]
%           of size L^2, where the first number corresponds to l
%           and the second to m.

if length(flm) ~= L^2
    error('Parameter flm has to contain L^2 coefficients.')
end

% step 1
% Sum of flm over l to obtain Fmm'
% 1st index is m', from -L+1 to L-1
% 2nd index is m, from -L+1 to L-1
fmm = zeros(2*L-1, 2*L-1);
m0i = L;

dl = zeros(2*L-1,2*L-1);
for l=0:L-1,
    dl = ssht_dl(dl, L, l, pi/2);
    
    lm0i = l^2+1+l;
    for m=-l:l,
        sign = 1i^(-m);
        for mp=-l:l, % m' in the equation (14)
            fmm(m0i+mp,m0i+m) = fmm(m0i+mp,m0i+m) +...
                sign*sqrt((2*l+1)/(4*pi))*...
                dl(m0i+mp,m0i+m)*dl(m0i+mp,m0i)*flm(lm0i+m) *...
                exp(1i*(mp)*pi/(2*L-1));
                % this phase shift is necessary to account for the desired
                % theta sample points (which are half a sample step size 
                % off 0)
        end
    end
end

% step 2
% 2d-IFFT of Fmm' over m' and m to obtain f (result)
% Because Fmm' does not have m' = 0 at the first position, we correct the 
% result using a phase shift.

phase = exp(2*pi*1i/(2*L-1)*(0:2*L-2)*L); % shifting by -(L-1) instead of 
                                          % +L would have the same effect

phasem = repmat(phase, L, 1);
phasemp = repmat(phase(1:L).', 1, 2*L-1); 

% 1st index is t, from 0 to 2L-2
% 2nd index is p, from 0 to 2L-2
f = ifft2(fmm);
% now drop additional theta values, scale, apply phase shifts
f = f(1:L,:).*(2*L-1)^2.*phasem.*phasemp; 