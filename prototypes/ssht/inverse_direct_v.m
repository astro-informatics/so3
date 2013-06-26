function [ f ] = inverse_direct_v( L, flm )
%INVERSE_DIRECT Inverse spherical harmonic transform using only matlab and
%a vectorised implementation.
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

c = cos(thetas);

for l = 0:L-1,
    ppos = legendre(l, c);
    if l > 0
        % Copy results for positive m to reuse them for negative m
        mneg = -l:-1;
        
        scale = repmat(((-1).^mneg.*factorial(l+mneg)./factorial(l-mneg))',1,size(ppos,2));
        pneg = scale.*ppos(end:-1:2,1:end);
        
        p = [pneg; ppos];
    else
        p = ppos;
    end
    % p is now a matrix containing all P(l,m,theta) for the current l
    % the 1st dimension corresponds to m, the 2nd on to the theta samples
    ms = -l:l;
    coeffs = flm(l^2+1:(l+1)^2);
    scale = sqrt((2*l+1)/(4*pi)*factorial(l-ms)./factorial(l+ms));
    scalem = repmat((coeffs.*scale)',1,size(p,2));
    e = exp(1i*(phis*(ms))); 
    % e contains all exponentials for all m and all phi
    % the 1st dimension corresponds to the phi samples, the 2nd to m
    f=f+(e*(scalem.*p)).';
end