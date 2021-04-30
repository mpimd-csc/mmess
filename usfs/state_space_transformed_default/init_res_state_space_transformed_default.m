function [RHS, res0, eqn, opts, oper] = ...
    init_res_state_space_transformed_default(eqn, opts, oper, RHS)
%% function [RHS, res0, eqn, opts, oper] = ...
%    init_res_state_space_transformed_default(eqn, opts, oper, RHS)
%
% This function returns the initial residual factor RHS and its
% associated norm res0.
%
% Input
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   oper            struct contains function handles for operation
%                   with A and E
%
%   RHS             right hand-side matrix
%
% Output
%   RHS             right hand-side matrix
%
%   res0            residuum norm of RHS
%
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   oper            struct contains function handles for operation
%                   with A and E

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Check input data.
if not(isnumeric(RHS)) || not(ismatrix(RHS))
    error( ...
        'MESS:error_arguments', ...
        'RHS has to ba a matrix');
end

%% Compute res0.
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
