function pflmn = so3_pack_flmn(flmn,L,N)

pflmn = zeros((N)*L*L, 1);

for el=0:L-1
    for em=-el:el
        n = min(el,N-1);
        for en=0:n
            ind_up = so3_elmn2ind(el, em, en, L, N, 'Order', 'NegativeFirst', 'Reality', false);
            ind_p  = so3_elmn2ind(el, em, en, L, N, 'Order', 'NegativeFirst', 'Reality', true);
            
            pflmn(ind_p) = flmn(ind_up);
        end
    end
end

