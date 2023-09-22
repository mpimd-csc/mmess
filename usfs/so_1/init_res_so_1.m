function [W, res0, eqn, opts, oper] = ...
    init_res_so_1(eqn, opts, oper, W, T)
%% function init_res initializes the low-rank residual W and res0
% function [ W, res0, eqn, opts, oper ] = ...
%      init_res_so_1( eqn, opts, oper, W, T)
%
% Call help mess_usfs_so_1 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns a matrix W and its associated residuum norm res0.
%
%   Input/Output:
%
%   eqn        structure containing data for G or B or C
%   opts       structure containing parameters for the algorithm
%   oper       struct contains function handles for operation with A and E
%   W          right hand side matrix
%   T          matrix such that the residual is W*T*W'
%              (optional, defaults to the identity)
%
%   Outputs:
%
%   W        matrix given by ADI to compute residuum
%   res0       initial residuum norm
%
% This function does not use other so1 functions.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input data
if (not(isnumeric(W))) || (not(ismatrix(W)))
    mess_err(opts, 'error_arguments', 'W has to ba a matrix');
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
