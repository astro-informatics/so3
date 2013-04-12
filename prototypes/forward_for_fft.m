function [ flm ] = forward_for_fft( L, f )
%FORWARD_FOR_FFT Forward spherical harmonic transform using matlab's fft
%and a factoring of rotations approach.
%.
%   This computes the forward spherical harmonic transform, synthesising 
%   the function using factoring of rotations. It also employs matlab's fft
%   to compute the DFT.
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
fm = zeros(2*L-1,2*L-1); % first index is m, second is t from -L+1 to L-1
fmm = zeros(2*L-1);
gmm = zeros(2*L-1, 2*L-1);

for t=0:length(thetas)-1,
    % step 1
    fm(:,L+t) = fftshift(fft(f(t+1,:)))./(2*L-1);
    % step 2
    if t < length(thetas)-1
        for m = -L+1:L-1,
            fm(L+m,L-t-1) = (-1)^m.*fm(L+m,L+t);
        end
    end
end

for m=-L+1:L-1,
    % step 3
    correction = exp(1i*pi*(0:2*L-2)/(2*L-1));
    fmm = fftshift(fft(ifftshift(fm(L+m,:)))./correction)./(2*L-1);
    % very bad hack... need to figure out the source of this discrepancy:
    fmm(1:L-1) = -fmm(1:L-1); 
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