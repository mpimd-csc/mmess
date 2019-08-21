function [ W, res0, eqn, opts, oper ] = init_res_dae_1_so( eqn, opts, oper, RHS)
%% function init_res initializes the low rank residual W and res0
% function [ W, res0, eqn, opts, oper ] = init_res_dae_1_so( eqn, opts, oper, RHS)
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
%   W          matrix given by ADI to compute residuum
%   res0       initial residuum norm
%
%   uses no other dae_1_so function
%% check data in eqn structure

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
%               2009-2019
%
if (not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_)))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
end
if (not(isfield(eqn,'isSym')))
    isSym = 0;
else
    isSym = eqn.isSym;
end
if not(isfield(eqn, 'nd'))    || not(isnumeric(eqn.nd))
    error('MESS:nd',...
    'Missing or Corrupted nd field detected in equation structure.');
end
if not(isfield(eqn,'type'))
  eqn.type='N';
  warning('MESS:equation_type',['Unable to determine type of equation.'...
    'Falling back to type ''N''']);
end

rowK = size(eqn.K_, 1);
nd = eqn.nd;
na = rowK - nd;

%% check input Paramters
if (not(isnumeric(RHS))) || (not(ismatrix(RHS)))
    error('MESS:error_arguments','RHS has to ba a matrix');
end
[rowG, colG] = size(RHS);


if(nd + na ~= rowG)
    error('MESS:error_arguments', 'number of rows of RHS'' has to coincide with dimension of M, D, K.');
end
%% compute low rank residual

if eqn.type == 'B'
    w = RHS(1 : nd, : ) - eqn.K_(1 : nd, nd + 1 : end) * ...
        (eqn.K_(nd + 1 : end, nd + 1 : end) \ (RHS(1 : nd, : )));
else
    if isSym
        w = RHS(1 : nd, : ) - eqn.K_(1 : nd, nd + 1 : end) * ...
            (eqn.K_(nd + 1 : end, nd + 1 : end) \ (RHS(1 : nd, : )));
    else
        w = RHS(1 : nd, : ) - eqn.K_(nd + 1 : end, 1 : nd)' * ...
            (eqn.K_(nd + 1 : end, nd + 1 : end)' \ (RHS(1 : nd, : )));
    end
end

W = [zeros(nd, colG); w];
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

