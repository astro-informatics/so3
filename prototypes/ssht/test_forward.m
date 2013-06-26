%   Tests prototype implementations of the forward transform.
%
%   This compares all matlab prototypes of the forward transform with the
%   ssht implementation and throws and error if the discrepancy exceeds
%   machine precision.

L = 8;
flm = random('unif', 0, 1, 1, L^2);
f = ssht_inverse(flm, L);
ssht = ssht_forward(f, L);
fine = 1;

result = forward_sov_for_dft(L, f);
err = max(abs(ssht(:) - result(:)))
if err > 10^-14
    disp(['Implementation forward_for_dft disagrees with ssht with errors up to ', num2str(err)]);
    fine = 0;
end

result = forward_sov_for_fft_phase(L, f);
err = max(abs(ssht(:) - result(:)))
if err > 10^-14
    disp(['Implementation forward_sov_for_fft_phase disagrees with ssht with errors up to ', num2str(err)]);
    fine = 0;
end

result = forward_sov_for_fft_space(L, f);
err = max(abs(ssht(:) - result(:)))
if err > 10^-14
    disp(['Implementation forward_sov_for_fft_space disagrees with ssht with errors up to ', num2str(err)]);
    fine = 0;
end

if fine == 1
    disp('All implementations agree with ssht up to machine precision.');
end

clear L
clear flm
clear f
clear ssht
clear fine
clear result
clear err