function [ flm ] = forward_sov_for_dft( L, f )
%FORWARD_FOR_DFT Forward spherical harmonic transform using only matlab and
%a factoring of rotations approach.
%.
%   This computes the forward spherical harmonic transform, synthesising 
%   the function using factoring of rotations.
%
%   L   ... the band limit (maximum l is L-1)
%   f   ... the sampled function values

[nt,np] = size(f);
if nt ~= L || np ~= 2*L-1
    error('Parameter f has to contain L * (2*L-1) coefficients.')
end

[thetas, phis] = ssht_sampling(L);

thetasExt = [- thetas(end-1:-1:1); thetas];

flm = zeros(L^2,1);
fm = zeros(2*L-1,1); % first index is m, second is t from -L+1 to L-1
fmm = zeros(2*L-1,1);
gmm = zeros(2*L-1, 2*L-1);

for m=-L+1:L-1,
    for t=0:length(thetas)-1,
        % step 1
        sum = 0;
        for p=1:length(phis),
            phi = phis(p);
            sum = sum + f(t+1,p)*exp(-1i*m*phi);
        end
        sum = sum/(2*L-1);
        fm(L+t) = sum;
        % step 2
        if t < length(thetas)-1
            fm(L-t-1) = (-1)^m*sum;
        end
    end

    % step 3
%     for mp=-L+1:L-1,
%         sum = 0;
%         for t=0:length(thetasExt)-1,
%             theta = thetasExt(t+1);
%             sum = sum + fm(t+1)*exp(-1i*mp*theta);
%         end
%         fmm(L+mp) = sum / (2*L-1);
%     end
    fm = ifftshift(fm);
    N = 2*L-1;
    for mp=1:N,
        if mp > L
            mpp = mp-N;
        else
            mpp = mp;
        end
            
        sum = 0;
        for t=1:N,
            sum = sum + fm(t)*exp(-2*pi*1i/N*(mp-1)*(t-1)); % here we can choose either mpp or mp
        end
        halfIntOffset = exp(-1i*(mpp-1)*pi/N); % here we have to use mpp
        fmm(mp) = sum*halfIntOffset / (2*L-1);
    end
    fmm = fftshift(fmm);
    
    % step 4
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
            sum = sum + fmm(L+mpp)*w;
        end
        gmm(L+m,L+mp) = sum * 2*pi;
    end
end

% step 5
dl = zeros(2*L-1,2*L-1);
lmi = 1;
for l=0:L-1,
    dl = ssht_dl(dl, L, l, pi/2);
    for m=-l:l,
        sum = 0;
        for mp=-L+1:L-1,
            sum = sum + dl(L+mp,L+m)*dl(L+mp,L)*gmm(L+m,L+mp);
        end
        flm(lmi) = sum * 1i^m * sqrt((2*l+1)/(4*pi));
        lmi = lmi + 1;
    end
end