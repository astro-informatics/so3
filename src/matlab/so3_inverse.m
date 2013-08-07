function [f] = so3_inverse(flmn, L, N, varargin)
% so3_inverse - Compute inverse Wigner transform
%
% Computes inverse Wigner transform.
%
% Default usage is given by
%
%   f = so3_inverse(flmn, L, N, <options>)
%
% where L and N are the harmonic and orientational band-limit,
% respectively, flmn is the vector of  harmonic coefficients and
% f is the sampled function values indexed by alpha, beta, gamma.
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%
%   'Order'   = { 'ZeroFirst'     [el-em blocks are stored in en-order
%                                  0, -1, 1, -2, 2, ... (default)],
%                 'NegativeFirst' [el-em blocks are stored in en-order
%                                  ..., -2, -1, 0, 1, 2, ...] }
%   'Storage' = { 'Padded'        [indices for el < en are zero (default)],
%                 'Compact'       [indices for el < en are omitted] }
%
% Note that an alternative interface is provided for the MW and MWSS
% methods since these contain samples on the poles (MW sampling
% contains a sample on the South pole, but not North pole; MWSS
% sampling contains a sample on both the North and South poles).
% Instead of requesting function values on the poles for all values of
% phi, the poles may be specified by single samples with their
% corresponding phi angle (values for other phi are then related to
% this sample by the rotation of a spin function in its tangent
% plane).  For the MW method this interface is called through the
% usage
%
%   [f, f_sp, phi_sp] = ssht_inverse(flm, L, <options>)
%
% where f does not contain samples on the South pole, f_sp is the
% South pole sample and phi_sp is its corresponding phi value.  For
% the MWSS method this interface is called through the usage
%
%   [f, f_sp, phi_sp, f_np, phi_np] = ssht_inverse(flm, L, <options>)
%
% where f_np is the North pole sample and phi_np is its corresponding
% phi value.
%
% Author: Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013  Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('flmn', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('N', @isnumeric);
p.addParamValue('Order', 'ZeroFirst', @ischar);
p.addParamValue('Storage', 'Padded', @ischar);
p.parse(flmn, L, N, varargin{:});
args = p.Results;

% Computing inverse transform.
[f] = ...
    so3_inverse_mex(flmn, L, N, args.Order, args.Storage);
