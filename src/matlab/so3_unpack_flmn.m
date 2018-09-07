function flmn = so3_unpack_flmn(pflmn,L,N)

flmn = zeros((2*N-1)*L*L, 1);

for el=0:L-1
    for em=-el:el
        n = min(el,N-1);
        for en=0:n
            ind_up_p = so3_elmn2ind(el, em, en,  L, N, 'Order', 'NegativeFirst', 'Reality', false);
            ind_up_m = so3_elmn2ind(el, -em, -en, L, N, 'Order', 'NegativeFirst', 'Reality', false);
            ind_p    = so3_elmn2ind(el, em, en,  L, N, 'Order', 'NegativeFirst', 'Reality', true);
            
            flmn(ind_up_p) = pflmn(ind_p);
            flmn(ind_up_m) = pflmn(ind_p)'*(-1)^(em-en);
        end
    end
end
