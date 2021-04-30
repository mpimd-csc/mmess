function [ RHS, res0, eqn, opts, oper ] = init_res_so_2( eqn, opts, oper, RHS)
%% function init_res initializes the low rank residual W and res0
% function [ RHS, res0, eqn, opts, oper ] = init_res_so_2( eqn, opts, oper, RHS)
%
% Call help mess_usfs_so_2 to see the description of the second order
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
%   RHS        right hand side matrix
%
%   Outputs:
%
%   RHS        matrix given by ADI to compute residuum
%   res0       initial residuum norm
%
% This function does not use other so3 functions.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input data
if (not(isnumeric(RHS))) || (not(ismatrix(RHS)))
    error('MESS:error_arguments','RHS has to ba a matrix');
end

%% compute res0
if isfield(opts, 'nm') && isfield(opts.nm, 'res0')
    res0 = opts.nm.res0;
else
    if opts.LDL_T
        if opts.norm == 2
            res0 = max(abs(eig(RHS' * RHS * diag(eqn.S_diag))));
        else
            res0 = norm(eig(RHS' * RHS * diag(eqn.S_diag)), 'fro');
        end
    else
        res0 = norm(RHS' * RHS, opts.norm);
    end
end

end

