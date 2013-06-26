%   Tests all prototypes and plots errors and computation times.

maxPower = 7;
time_1d = zeros(maxPower,1);
time_md_col = zeros(maxPower,1);
time_md_row = zeros(maxPower,1);
time_no_array = zeros(maxPower,1);

so3_disagreement = zeros(maxPower,1);

for Lpower = 1:maxPower,
    L = 2^Lpower;
    M = L;
    N = L;
    
    flmn = random('unif', 0, 1, [(2*N-1)*(3*L^2-N*(N-1))/3, 1]);
    
    tic;
    f = so3_direct_inverse_1d(L, M, N, flmn);
    time_1d(Lpower) = toc;
    
    tic;
    f = so3_direct_inverse_md(L, M, N, flmn);
    time_md_col(Lpower) = toc;
    
    tic;
    f = so3_direct_inverse(L, M, N, flmn);
    time_md_row(Lpower) = toc;
    
    tic;
    f = so3_direct_inverse_no_array(L, M, N, flmn);
    time_no_array(Lpower) = toc;
end

x = 1:maxPower;

figure
semilogy(x, 2.^(4:4:4*maxPower).'./1000, '-k',...
         x, time_1d,     '-r',...
         x, time_md_col,   '-b',...
         x, time_md_row,  '-g',...
         x, time_md_row,  '-m');
title('Computation time of inverse transform');
legend('O(L^4)', '1d-access', '3d-access by columns', '3d-access by "rows"', 'no array access', 'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

clear L M N maxPower
clear flmn f result
clear time_1d time_md_col time_md_row time_no_array