function [W, res0, eqn, opts, oper] = ...
    init_res_default_iter(eqn, opts, oper, W, T)
%% Function init_res_default_iter initializes the low-rank residual W and res0
% function [ W, res0, eqn, opts, oper ] = ...
%     init_res_default_iter( eqn, opts, oper, W, T)
%
% This function returns the initial residual factor W and its
% associated norm res0.
%
%   Input/Output:
%
%   eqn        structure containing data for G or B or C
%   opts       structure containing parameters for the algorithm
%   oper       structure contains function handles for operation with A and E
%   W          right hand side matrix
%   T          matrix such that the residual is W*T*W'
%              (optional, defaults to the identity)
%
%   Outputs:
%
%   W        matrix given by ADI to compute residuum
%   res0       initial residuum norm
%
% This function does not use other default functions.
%

% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
%% Check input data
if (not(isnumeric(W))) || (not(ismatrix(W)))
    mess_err(opts, 'error_arguments', 'W has to ba a matrix');
end

%% Compute res0
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
