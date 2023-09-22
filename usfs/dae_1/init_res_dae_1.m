function [W, res0, eqn, opts, oper] = ...
    init_res_dae_1(eqn, opts, oper, W, T)
%% function init_res initializes the low-rank residual W and res0
% function [ W, res0, eqn, opts, oper ] = ...
%          init_res_dae_1( eqn, opts, oper, W, T)
%
%   Input/Output:
%
%   eqn        structure containing data for G or B or C
%   opts       structure containing parameters for the algorithm
%   oper       struct contains function handles for operation with A and E
%   W          right hand side matrix
%   T          such that residual is W*T*W'
%
%   Outputs:
%
%   W        matrix given by ADI to compute residuum
%   res0       initial residuum norm
%
%   uses no other dae_1 function

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check data
if not(isfield(opts, 'LDL_T'))
    opts.LDL_T = false;
end
if opts.LDL_T && not(isfield(eqn, 'T'))
    mess_err(opts, 'control_data', 'T required in LDL_T mode');
end

if (not(isnumeric(W))) || (not(ismatrix(W)))
    mess_err(opts, 'error_arguments', 'W has to ba a matrix');
end
if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'control_data', ...
             'Missing or corrupted manifold_dim field detected in eqn structure.');
end
if not(eqn.manifold_dim == size(W, 1))
    mess_err(opts, 'error_arguments', ...
             'number of rows of A_ differs with number of rows of W');
end

%% compute res0
if not(exist('T', 'var')) && opts.LDL_T
    % this means we only use init_res for potential projection
    return
end
if isfield(opts, 'nm') && isfield(opts.nm, 'res0')
    res0 = opts.nm.res0;
else
    if opts.LDL_T
        if opts.norm == 2
            res0 = max(abs(eig(W' * W * T)));
        else
            res0 = norm(eig(W' * W * T), 'fro');
        end
    else
        res0 = norm(W' * W, opts.norm);
    end
end

end
