function [el, m] = ssht_ind2elm(ind)
% ssht_ind2elm - Convert vector index to harmonic indices
%
% Convert ind index to access vector of harmonic coefficients (folllowing
% the Matlab convention, ind is index from 1) to (el,m) spherical harmonic
% indices. 
%
% Default usage is given by
%
%   [el, m] = ssht_ind2elmn(ind)
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SSHT package to perform spin spherical harmonic transforms
% Copyright (C) 2011  Jason McEwen
% See LICENSE.txt for license details

el = floor(sqrt(ind-1));
m = ind - 1 - el*el - el;
