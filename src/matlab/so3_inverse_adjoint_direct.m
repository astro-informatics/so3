function [flmn] = so3_inverse_adjoint_direct(f, L, N, varargin)
% so3_forward - Compute adjoint inverse Wigner transform
%
% Computes adjoint inverse Wigner transform.
%
% Default usage is given by
%
%   flmn = so3_inverse_adjoint_direct(f, L, N, <options>)
%
% where L and N are the harmonic and orientational band-limit,
% respectively, flmn is the vector of  harmonic coefficients and
% f is the sampled function values indexed by gamma, beta, alpha.
%
% Options consist of parameter type and value pairs.  Valid options
% include:
%
%   'L0' = Non-negative integer (default 0). Indicates that only flmn
%          with l >= L0 are non-zero. Currently there are no savings
%          in storage but there are some in performance.
%
%   'Sampling' = { 'MW'           [McEwen & Wiaux sampling (default)],
%                  'MWSS'         [McEwen & Wiaux symmetric sampling] }
%
%   'Order' = { 'ZeroFirst'       [el-em blocks are stored in en-order
%                                  0, -1, 1, -2, 2, ... (default)],
%               'NegativeFirst'   [el-em blocks are stored in en-order
%                                  ..., -2, -1, 0, 1, 2, ...] }
%
%   'Storage' = { 'Padded'        [indices for el < en are zero (default)],
%                 'Compact'       [indices for el < en are omitted] }
%
%   'NMode' = { 'All'             [all flmn are non-zero (default)],
%               'L'               [only flmn for |m| = L are non-zero],
%               'Even'            [only flmn for even n are non-zero],
%               'Odd'             [only flmn for odd n are non-zero],
%               'Maximum'         [only flmn for |n| = N-1 are non-zero] }
%
%   'WignerMethod' = { 'Trapani'  [Trapani-recursion is used to compute the
%                                  the Wigner functions (default)],
%                      'Risbo'    [Risbo-recursion is used to compute the
%                                  the Wigner functions] }
%
%  'Reality' = { false            [do not assume f real (default)],
%                true             [assume f real (improves performance)] }
%
% Authors: Martin Büttner (m.buettner.d@gmail.com)
%          Jason McEwen (www.jasonmcewen.org)

% SO3 package to perform Wigner transforms
% Copyright (C) 2013 Martin Büttner and Jason McEwen
% See LICENSE.txt for license details

% Parse arguments.
p = inputParser;
p.addRequired('f', @isnumeric);
p.addRequired('L', @isnumeric);
p.addRequired('N', @isnumeric);
p.addParamValue('L0', 0, @isnumeric);
p.addParamValue('Sampling', 'MW', @ischar);
p.addParamValue('Order', 'ZeroFirst', @ischar);
p.addParamValue('Storage', 'Padded', @ischar);
p.addParamValue('NMode', 'All', @ischar);
p.addParamValue('WignerMethod', 'Trapani', @ischar);
p.addParamValue('Reality', false, @islogical);
p.parse(f, L, N, varargin{:});
args = p.Results;

% Computing adjoint inverse transform.
[flmn] = ...
    so3_inverse_adjoint_direct_mex(f, args.L0, L, N, args.Order, args.Storage, args.NMode, args.WignerMethod, args.Reality, args.Sampling);
