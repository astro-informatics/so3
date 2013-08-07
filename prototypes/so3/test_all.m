%   Tests all prototypes and plots errors and computation times.

maxPower = 5;
so3_direct_error = zeros(maxPower,1);
so3_direct_time_forward = zeros(maxPower,1);
so3_direct_time_inverse = zeros(maxPower,1);

so3_via_ssht_error = zeros(maxPower,1);
so3_via_ssht_time_forward = zeros(maxPower,1);
so3_via_ssht_time_inverse = zeros(maxPower,1);

so3_c_error = zeros(maxPower,1);
so3_c_time_forward = zeros(maxPower,1);
so3_c_time_inverse = zeros(maxPower,1);

so3_disagreement = zeros(maxPower,1);

for Lpower = 1:maxPower,
    L = 2^Lpower;
    M = L;
    N = L;
    
    flmn = random('unif', 0, 1, [(2*N-1)*(3*L^2-N*(N-1))/3, 1]);
    
    tic;
    f_direct = so3_direct_inverse(L, M, N, flmn);
    so3_direct_time_inverse(Lpower) = toc;
    tic;
    result = so3_direct_forward(L, M, N, f_direct);
    so3_direct_time_forward(Lpower) = toc;
    so3_direct_error(Lpower) = max(abs(flmn(:) - result(:)));
    
    tic;
    f_via_ssht = so3_via_ssht_inverse(L, N, flmn);
    so3_via_ssht_time_inverse(Lpower) = toc;
    tic;
    result = so3_via_ssht_forward(L, N, f_via_ssht);
    so3_via_ssht_time_forward(Lpower) = toc;
    so3_via_ssht_error(Lpower) = max(abs(flmn(:) - result(:)));
    
    so3_disagreement(Lpower) = max(abs(f_direct(:) - f_via_ssht(:)));
    
    tic;
    f_c = so3_inverse(flmn, L, N, 'Order', 'NegativeFirst', 'Storage', 'Compact');
    so3_c_time_inverse(Lpower) = toc;
    tic;
    result = so3_forward(f_c, L, N, 'Order', 'NegativeFirst', 'Storage', 'Compact');
    so3_c_time_forward(Lpower) = toc;
    so3_c_error(Lpower) = max(abs(flmn(:) - result(:)));
end

x = 1:maxPower;
figure
semilogy(x, so3_direct_error,   '-b',...
         x, so3_via_ssht_error, '-r',...
         x, so3_disagreement,   '-m',...
         x, so3_c_error,        '-g');
title('Errors in accuracy.');
legend('Round-trip error - direct implementation',...
       'Round-trip error - implementation via SSHT',...
       'Disagreement in f between implementations',...
       'Round-trip error - C implementation via SSHT',...
       'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

figure
semilogy(x, 2.^(4:4:4*maxPower).'./1000, '-k',...
         x, so3_direct_time_forward,     '-b',...
         x, so3_via_ssht_time_forward,   '-r',...
         x, so3_c_time_forward,          '-g');
title('Computation time of forward transform');
legend('O(L^4)', ...
       'Prototype - direct implementation', ...
       'Prototype - implementation via SSHT', ...
       'C implementation - via SSHT', ...
       'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

figure
semilogy(x, 2.^(4:4:4*maxPower).'./1000, '-k',...
         x, so3_direct_time_inverse,     '-b',...
         x, so3_via_ssht_time_inverse,   '-r',...
         x, so3_c_time_inverse,   '-g');
title('Computation time of inverse transform');
legend('O(L^4)', ...
       'Prototype - direct implementation', ...
       'Prototype - implementation via SSHT', ...
       'C implementation - via SSHT', ...
       'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

clear L M N maxPower
clear flmn f_direct f_via_ssht f_c result
clear so3_direct_error so3_direct_time_forward so3_direct_time_inverse
clear so3_via_ssht_error so3_via_ssht_time_forward so3_via_ssht_time_inverse
clear so3_c_error so3_c_time_forward so3_c_time_inverse
clear so3_disagreement