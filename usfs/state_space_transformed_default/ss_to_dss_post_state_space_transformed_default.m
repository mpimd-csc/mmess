function [eqn, opts, oper] = ...
    ss_to_dss_post_state_space_transformed_default(eqn, opts, oper)
%% function [eqn, opts, oper] = ...
%     ss_to_dss_post_state_space_transformed_default(eqn, opts, oper)
%
% function post finalizes data and/or functions
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

assert(isfield(eqn, 'Ecount'), ...
    'MESS:error_arguments', ...
    'field eqn.Scount is not defined.');

if eqn.Ecount > 1
    eqn.Ecount = eqn.Ecount - 1;
else
    eqn = rmfield(eqn, 'Ecount');
    eqn = rmfield(eqn, 'EL');
    eqn = rmfield(eqn, 'EU');
end
