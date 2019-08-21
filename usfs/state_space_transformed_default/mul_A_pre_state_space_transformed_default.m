function [eqn, opts, oper] = ...
    mul_A_pre_state_space_transformed_default(eqn, opts, oper)
%% function [eqn, opts, oper] = ...
%     mul_A_pre_state_space_transformed_default(eqn, opts, oper)
%
% function pre initializes data and/or functions
%
% Input
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   oper            struct contains function handles for operation
%                   with A and E
%
% Output
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
%               2009-2019
%

if isfield(eqn, 'EL') && isfield(eqn, 'EU')
    if isfield(eqn, 'Ecount')
        eqn.Ecount = eqn.Ecount + 1;
    else
        eqn.Ecount = 2;
    end
else
    [eqn.EL, eqn.EU] = lu(eqn.E_);
    eqn.Ecount       = 1;
end
