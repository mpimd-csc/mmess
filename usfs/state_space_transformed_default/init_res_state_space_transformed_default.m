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
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others
%               2009-2020
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
