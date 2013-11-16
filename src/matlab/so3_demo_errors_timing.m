% so3_demo_errors_timing - Create error and timing plots.
%
% Tests the C implementation and plots errors and computation times.
%
% Default usage is given by
%
%   so3_demo_errors_timing
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013 Martin Büttner and Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

maxPower = 7;

so3_NL_error = zeros(maxPower,1);
so3_NL_time_forward = zeros(maxPower,1);
so3_NL_time_inverse = zeros(maxPower,1);

so3_N3_error = zeros(maxPower,1);
so3_N3_time_forward = zeros(maxPower,1);
so3_N3_time_inverse = zeros(maxPower,1);

for Lpower = 1:maxPower,
    L = 2^Lpower;

    N = min(L,3);
    flmn = random('unif', 0, 1, [(2*N-1)*(3*L^2-N*(N-1))/3, 1]);

    tic;
    f = so3_inverse(flmn, L, N, 'Storage', 'Compact');
    so3_N3_time_inverse(Lpower) = toc;
    tic;
    result = so3_forward(f, L, N, 'Storage', 'Compact');
    so3_N3_time_forward(Lpower) = toc;
    so3_N3_error(Lpower) = max(abs(flmn(:) - result(:)));

    N = L;
    flmn = random('unif', 0, 1, [(2*N-1)*(3*L^2-N*(N-1))/3, 1]);

    tic;
    f = so3_inverse(flmn, L, N, 'Storage', 'Compact');
    so3_NL_time_inverse(Lpower) = toc;
    tic;
    result = so3_forward(f, L, N, 'Storage', 'Compact');
    so3_NL_time_forward(Lpower) = toc;
    so3_NL_error(Lpower) = max(abs(flmn(:) - result(:)));
end

x = 1:maxPower;
figure
semilogy(x, so3_N3_error, '-g',...
         x, so3_NL_error, '-b');
title('Errors in accuracy.');
legend('C implementation - N = 3',...
       'C implementation - N = L',...
       'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

figure
semilogy(x, 2.^(4:4:4*maxPower).'./10^6.5, '-k',...
         x, 2.^(3:3:3*maxPower).'./10^7.5, '--k',...
         x, so3_N3_time_forward,         '-g',...
         x, so3_NL_time_forward,         '-b');
title('Computation time of forward transform');
legend('O(L^4)', ...
       'O(L^3)', ...
       'C implementation - N = 3', ...
       'C implementation - N = L', ...
       'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

figure
semilogy(x, 2.^(4:4:4*maxPower).'./10^6.5, '-k',...
         x, 2.^(3:3:3*maxPower).'./10^7.5, '--k',...
         x, so3_N3_time_inverse,         '-g',...
         x, so3_NL_time_inverse,         '-b');
title('Computation time of inverse transform');
legend('O(L^4)', ...
       'O(L^3)', ...
       'C implementation - N = 3', ...
       'C implementation - N = L', ...
       'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

clear L N maxPower Lpower
clear flmn f result
clear x
clear so3_NL_error so3_NL_time_forward so3_NL_time_inverse
clear so3_N3_error so3_N3_time_forward so3_N3_time_inverse
