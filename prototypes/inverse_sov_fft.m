function [ f ] = inverse_sov_fft( L, flm )
%INVERSE_SOV_FFT Inverse spherical harmonic transform using matlab's ifft
%and a separation of variables approach.
%.
%   This computes the inverse spherical harmonic transform, synthesising 
%   the function using separation of variables. It also employs matlab's
%   ifft to compute the DFT.
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

[thetas, phis] = ssht_sampling(L);

f = zeros(length(thetas), length(phis));
fm = zeros(length(thetas), 2*L-1);
m0i = L;

for j=1:length(thetas),
    theta = thetas(j);
    dln = zeros(L,1);
    dlnprev = zeros(L,1);
    for l=0:L-1,
        tmp = ssht_dln(dln, dlnprev, L, l, 0, theta);
        dlnprev = dln;
        dln = tmp;
        
        lm0i = l^2+1+l;
        sign = (-1)^l;
        for m=-l:-1,
            fm(j, m0i+m) = fm(j, m0i+m) + flm(lm0i+m)*sqrt((2*l+1)/(4*pi))*sign*dln(abs(m)+1);
            sign = -sign;
        end
        for m=0:l,
            fm(j, m0i+m) = fm(j, m0i+m) + flm(lm0i+m)*sqrt((2*l+1)/(4*pi))*dln(abs(m)+1);
        end
    end
end

for j=1:length(thetas),
    theta = thetas(j);
    f(j,:) = ifft(fm(j,:)).*(2*L-1).*exp(-1i.*(L-1).*phis.');
end
