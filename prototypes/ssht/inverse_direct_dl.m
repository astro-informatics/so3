function [ f ] = inverse_direct_dl( L, flm )
%INVERSE_DIRECT Inverse spherical harmonic transform using ssht_dl.
%   This computes the inverse spherical harmonic transform, simply
%   synthesising the function by summing over all spherical harmonics.
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

for j=1:length(thetas),
    theta = thetas(j);
    dl = zeros(2*L-1,2*L-1);
    for l = 0:L-1,
        dl = ssht_dl(dl, L, l, theta);
        for k=1:length(phis),
            phi = phis(k);
            lmindex = l^2+1;
            for m = -l:l,
                f(j,k) = f(j,k) + flm(lmindex)*sqrt((2*l+1)/(4*pi))*dl(m+L, L)*exp(1i*m*phi);
                    lmindex = lmindex + 1;
            end
        end
    end
end