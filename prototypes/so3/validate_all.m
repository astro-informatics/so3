%   Tests all prototypes and plots errors and computation times.

maxPower = 6;

so3_inv_disagreement_pdirect_pssht = zeros(maxPower,1);
so3_inv_disagreement_pdirect_c = zeros(maxPower,1);
so3_inv_disagreement_pssht_c = zeros(maxPower,1);
so3_forw_disagreement_pdirect_pssht = zeros(maxPower,1);
so3_forw_disagreement_pdirect_c = zeros(maxPower,1);
so3_forw_disagreement_pssht_c = zeros(maxPower,1);

for Lpower = 1:maxPower,
    L = 2^Lpower;
    M = L;
    N = min(3,L);
    
    flmn = random('unif', 0, 1, [(2*N-1)*(3*L^2-N*(N-1))/3, 1]);
    
    f_direct = so3_direct_inverse(L, M, N, flmn);
    f_via_ssht = so3_via_ssht_inverse(L, N, flmn);
    f_c = permute(so3_inverse(flmn, L, N, 'Order', 'NegativeFirst', 'Storage', 'Compact'), [3 2 1]);
    
    so3_inv_disagreement_pdirect_pssht(Lpower) = max(abs(f_direct(:) - f_via_ssht(:)));
    so3_inv_disagreement_pdirect_c(Lpower) = max(abs(f_direct(:) - f_c(:)));
    so3_inv_disagreement_pssht_c(Lpower) = max(abs(f_via_ssht(:) - f_c(:)));
    
    flmn_direct = so3_direct_forward(L, M, N, f_c);
    flmn_via_ssht = so3_via_ssht_forward(L, N, f_c);
    flmn_c = so3_forward(permute(f_c, [3 2 1]), L, N, 'Order', 'NegativeFirst', 'Storage', 'Compact');
    
    so3_forw_disagreement_pdirect_pssht(Lpower) = max(abs(flmn_direct(:) - flmn_via_ssht(:)));
    so3_forw_disagreement_pdirect_c(Lpower) = max(abs(flmn_direct(:) - flmn_c(:)));
    so3_forw_disagreement_pssht_c(Lpower) = max(abs(flmn_via_ssht(:) - flmn_c(:)));
end

x = 1:maxPower;
figure
semilogy(x, so3_inv_disagreement_pdirect_pssht, '-b',...
         x, so3_inv_disagreement_pdirect_c,     '-r',...
         x, so3_inv_disagreement_pssht_c,       '-g');
title('Disagreement between inverse transform implementations.');
legend('Prototype (direct) vs prototype (via SSHT)',...
       'Prototype (direct) vs C implementation',...
       'Prototype (via SSHT) vs C implementation',...
       'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

figure
semilogy(x, so3_forw_disagreement_pdirect_pssht, '-b',...
         x, so3_forw_disagreement_pdirect_c,     '-r',...
         x, so3_forw_disagreement_pssht_c,       '-g');
title('Disagreement between forward transform implementations.');
legend('Prototype (direct) vs prototype (via SSHT)',...
       'Prototype (direct) vs C implementation',...
       'Prototype (via SSHT) vs C implementation',...
       'Location', 'NorthWest');
set(gca, 'XTick', (1:maxPower));
set(gca, 'XTickLabel', 2.^(1:maxPower));

clear L M N maxPower Lpowers
clear flmn f_direct f_via_ssht f_c result flmn_direct flmn_via_ssht flmn_c
clear x
clear so3_inv_disagreement_pdirect_pssht so3_forw_disagreement_pdirect_pssht
clear so3_inv_disagreement_pdirect_c so3_forw_disagreement_pdirect_c
clear so3_inv_disagreement_pssht_c so3_forw_disagreement_pssht_c