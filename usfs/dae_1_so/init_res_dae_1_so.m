function [ RHS, res0, eqn, opts, oper ] = init_res_dae_1_so( eqn, opts, oper, RHS)
%% function init_res initializes the low rank residual W and res0
% function [ RHS, res0, eqn, opts, oper ] = init_res_dae_1_so( eqn, opts, oper, RHS)
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
%   RHS        restriction of RHS to the hidden manifold
%   res0       initial residuum norm
%
%   uses no other dae_1_so function

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check data

if not(isfield(opts, 'LDL_T'))
    opts.LDL_T = 0;
end
if opts.LDL_T && not(isfield(eqn, 'S_diag'))
    error('MESS:control_data', 'eqn.S_diag required in LDL_T mode');
end

if (not(isnumeric(RHS))) || (not(ismatrix(RHS)))
    error('MESS:error_arguments','RHS has to ba a matrix');
end

if not(isfield(eqn, 'nd')) || not(isnumeric(eqn.nd))
    error('MESS:control_data',...
    'Missing or corrupted nd field detected in eqn structure.');
end

if (2*eqn.nd ~= size(RHS, 1))
    error('MESS:error_arguments', ...
        'number of rows of A_ differs with number of rows of RHS');
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

