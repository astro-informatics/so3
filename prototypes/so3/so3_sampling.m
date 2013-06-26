function [ alphas, betas, gammas ] = so3_sampling( L, M, N )
%SO3_SAMPLING Compute sample positions
% 
%  Computes sample positions on sphere for various sampling methods.
% 
%  Default usage is given by
% 
%    [alphas, betas, gammas] = so3_sampling(L, M, N)
% 
%  where L, M and N are the harmonic band-limits for alpha, beta and gamma,
%  respectively, and alphas, betas and gammas specify sample positions. 

alphas = 2*(0:2*M-2).'*pi/(2*M-1);
betas  = (2*(0:L-1)+1).'*pi/(2*L-1);
gammas = 2*(0:2*N-2).'*pi/(2*N-1);

end

