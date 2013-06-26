function [ f ] = inverse_direct_dln( L, flm )
%INVERSE_DIRECT Inverse spherical harmonic transform using ssht_dln.
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
    dln = zeros(L,1);
    dlnprev = zeros(L,1);
    for l=0:L-1,
        tmp = ssht_dln(dln, dlnprev, L, l, 0, theta);
        dlnprev = dln;
        dln = tmp;
        for k=1:length(phis),
            phi = phis(k);
            lmindex = l^2+1;            
            for m = -l:l,
                sign = 1;
                if m < 0
                    sign = (-1)^m;
                end
                f(j,k) = f(j,k) + flm(lmindex)*sqrt((2*l+1)/(4*pi))*sign*dln(abs(m)+1)*exp(1i*m*phi);
                lmindex = lmindex + 1;
            end
        end
    end
end