function [ind] = so3_elmn2ind(el, em, en, L, N, varargin)
% so3_elmn2ind - Convert harmonic indices to vector index
%
% Convert (el,em,en) spherical harmonic indices to ind index to access vector
% of harmonic coefficients (following the Matlab convention, ind is index
% from 1).
%
% Default usage is given by
%
%   [ind] = so3_elmn2ind(el, em, en, L, N, <options>)
%
% Where el, em and en are the harmonic, azimuthal and orientational index,
% respectively and L and N are the harmonic and orientational band-limit,
% respectively.
%
% Options consist of parameter name and value pairs. Valid options include:
%
%   'Order'   = { 'ZeroFirst'     [el-em blocks are stored in en-order
%                                  0, -1, 1, -2, 2, ... (default)],
%                 'NegativeFirst' [el-em blocks are stored in en-order
%                                  ..., -2, -1, 0, 1, 2, ...] }
%   'Storage' = { 'Padded'        [indices for el < en are zero (default)],
%                 'Compact'       [indíces for el < en are omitted] }
%
%  'Reality' = { false            [do not assume f real (default)],
%                true             [assume f real; in this case, the Order
%                                  parameter will be ignored] }

% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013 Martin Büttner and Jason McEwen
% See LICENSE.txt for license details

% Parse arguments
p = inputParser;
p.addRequired('el', @isnumeric);
p.addRequired('em', @isnumeric);
p.addRequired('en', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('N', @isnumeric);
p.addParamValue('Order', 'ZeroFirst', @ischar);
p.addParamValue('Storage', 'Padded', @ischar);
p.addParamValue('Reality', false, @islogical);
p.parse(el, em, en, L, N, varargin{:});
args = p.Results;

% Compute index
[ind] = so3_elmn2ind_mex(el, em, en, L, N, args.Order, args.Storage, args.Reality);
