function [ eqn, opts, oper ] = eval_matrix_functions_default( eqn, opts, oper, t )
%% function eval_matrix_functions_default updates the matrices in eqn

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
if eqn.LTV
    
    %%
    if eqn.haveE
        eqn.E_ = eqn.E_time(t);
        eqn.A_ = eqn.A_time(t) + eqn.dt_E_time(t);
    else
        eqn.A_ = eqn.A_time(t);
    end
    eqn.B = eqn.B_time(t);
    eqn.C = eqn.C_time(t);
end
end