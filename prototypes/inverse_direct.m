function [ f ] = inverse_direct( L, flm )
%INVERSE_DIRECT Inverse spherical harmonic transform using naive approach
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
    for k=1:length(phis),
        phi = phis(k);
        sum = 0;
        for l = 0:L-1,
            ppos = legendre(l,cos(theta));
            for m = -l:l,
                if m >= 0
                    p = ppos(m+1);
                else
                    p = (-1)^m*factorial(l+m)/factorial(l-m)*ppos(-m+1);
                end
                sum = sum + flm(l^2 + l + 1 + m)*...
                    sqrt(...
                      (2*l+1)*factorial(l-m)/...
                      ((4*pi)*factorial(l+m))...
                    )*...
                    p*exp(1i*m*phi);
            end
        end
        f(j,k) = sum;
    end
end