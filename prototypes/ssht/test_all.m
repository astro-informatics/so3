%   Tests all prototypes and plots errors and computation times.

maxPower = 4;
ssht_error = zeros(maxPower,1);
ssht_time_forward = zeros(maxPower,1);
ssht_time_inverse = zeros(maxPower,1);

sov_for_phase_error = zeros(maxPower,1);
sov_for_phase_time_forward = zeros(maxPower,1);
sov_for_phase_time_inverse = zeros(maxPower,1);

sov_for_space_error = zeros(maxPower,1);
sov_for_space_time_forward = zeros(maxPower,1);
sov_for_space_time_inverse = zeros(maxPower,1);

for Lpower = 1:maxPower,
    L = 2^Lpower;
    flm = random('unif', 0, 1, 1, L^2);
    
    tic;
    f = ssht_inverse(flm, L);
    ssht_time_inverse(Lpower) = toc;
    tic;
    result = ssht_forward(f, L);
    ssht_time_forward(Lpower) = toc;
    ssht_error(Lpower) = max(abs(flm(:) - result(:)));
    
    tic;
    f = inverse_sov_for_fft_phase(L, flm);
    sov_for_phase_time_inverse(Lpower) = toc;
    tic;
    result = forward_sov_for_fft_phase(L, f);
    sov_for_phase_time_forward(Lpower) = toc;
    sov_for_phase_error(Lpower) = max(abs(flm(:) - result(:)));
    
    tic;
    f = inverse_sov_for_fft_space(L, flm);
    sov_for_space_time_inverse(Lpower) = toc;
    tic;
    result = forward_sov_for_dft(L, f);
    sov_for_space_time_forward(Lpower) = toc;
    sov_for_space_error(Lpower) = max(abs(flm(:) - result(:)));    
end

x = 1:maxPower;
figure
semilogy(x, ssht_error,          '-b',...
         x, sov_for_phase_error, '-r',...
         x, sov_for_space_error, '-g');
title('Errors in accuracy of round-trip transform');
legend('SSHT', 'Prototype - Phase', 'Prototype - Space', 'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));
 
figure
semilogy(x, 2.^(3:3:3*maxPower).'./1000, '-k',...
         x, ssht_time_forward,           '-b',...
         x, sov_for_phase_time_forward,  '-r',...
         x, sov_for_space_time_forward,  '-g');
title('Computation time of forward transform');
legend('O(L^3)', 'SSHT', 'Prototype - Phase', 'Prototype - Space', 'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

figure
semilogy(x, 2.^(3:3:3*maxPower).'./1000, '-k',...
         x, ssht_time_inverse,           '-b',...
         x, sov_for_phase_time_inverse,  '-r',...
         x, sov_for_space_time_inverse,  '-g');
title('Computation time of inverse transform');
legend('O(L^3)', 'SSHT', 'Prototype - Phase', 'Prototype - Space', 'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

clear L maxPower
clear flm f result
clear ssht_error ssht_time_forward ssht_time_inverse
clear sov_for_phase_error sov_for_phase_time_forward sov_for_phase_time_inverse
clear sov_for_space_error sov_for_space_time_forward sov_for_space_time_inverse