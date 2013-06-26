function [ f ] = inverse_sov_for_fft_space( L, flm )
%INVERSE_FOR_FFT Inverse spherical harmonic transform using matlab's ifft
%and a factoring of rotations approach.
%.
%   This computes the inverse spherical harmonic transform, synthesising 
%   the function using factoring of rotations (as outlined in equations
%   (12) to (14) of the ssht paper). It also employs matlab's ifft to 
%   compute the DFT.
%   The discrepancy between the required transform and matlab's definition
%   is treated by a spatial shift of the frequencies in m and m'. As the
%   function values are sampled at points in theta which are half a sample
%   step off zero, there is also a phase shift in the m' to account for
%   this spatial shift.
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
% But first, shift the 0-frequency component (m = m' = 0) to the top left
% corner of Fmm' (as matlab expects it)
f = ifft2(ifftshift(fmm));

% drop additional theta values and scale up
f = f(1:L,:).*(2*L-1)^2; 