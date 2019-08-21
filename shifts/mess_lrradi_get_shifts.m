function [eqn, opts, oper, nShifts] = ...
    mess_lrradi_get_shifts(eqn, opts, oper, W, Z, Y)
% Compute the next batch of shifts for the RADI method.
% 
% Input:
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation 
%                       with A and E
%
%   W                   the current residual matrix in RADI
%
%   Z                   the current Z matrix in RADI
%
%   Y                   the small square factor in the RADI solution
%
% opts.shifts.method 	possible  values:
%                           'precomputed'
%                           'penzl'
%                           'projection'
%                           'gen-ham-opti'
%
% If opts.shifts.method == 'heur', then:
% - calls mess_para to generate the shifts.
%
% If opts.shifts.method == 'projection', then:
% - calls mess_get_projection_shifts to generate the shifts.
%
% If opts.shifts.method == 'gen-ham-opti', then:
% - calls mess_lrradi_get_shifts_hamOpti_generalized to generate the shifts.
%

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
%               2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2018
%

%% Check input
switch opts.shifts.method
    case 'precomputed'
        % Just return the same array of shifts as the user provided.
        nShifts = length(opts.shifts.p);

    case 'heur' 
        % Use MESS routines for heuristic penzl shifts.
        p             = mess_para(eqn, opts, oper);
        nShifts       = length(p);
        opts.shifts.p = p;
    
    case 'projection' 
        % Use MESS routines for projection shifts.
        if isempty(Z)
            Z = W; 
        end
        
        nShifts = 0;
        i = 1;
        while nShifts==0
            nZ = size(Z, 2);
            if not(opts.radi.compute_sol_facpart)
                % Relevant parts of LRF.
                maxcolZ = opts.shifts.history;
                ind     = max(nZ - maxcolZ, 0)+1:nZ;
            else
                ind = 1:nZ;
            end
            
            if eqn.type == 'T'
                p = mess_projection_shifts( eqn, opts, oper, Z(:,ind), ...
                    oper.mul_A(eqn, opts, 'T', Z(:,ind), 'N') ...
                    + eqn.V * (eqn.U' * Z(:,ind)), []);
            else
                p = mess_projection_shifts( eqn, opts, oper, Z(:,ind), ...
                    oper.mul_A(eqn, opts, 'N', Z(:,ind), 'N') ...
                    + eqn.U * (eqn.V' * Z(:,ind)), []);
            end
            
            nShifts       = length(p);
            opts.shifts.p = p;
            if isempty(p)
                if  (i < 5)
                    warning('MESS:mess_para',...
                        ['Could not compute initial projection shifts. ',...
                        'Going to retry with random RHS.']);
                    Z = rand(size(Z));
                else
                    error('MESS:mess_para',...
                        'Could not compute initial projection shifts.');
                end
            end
            i = i + 1;
        end
    case 'gen-ham-opti'
        % Residual hamiltonian shifts... recommended+default!
        [eqn, opts, oper, nShifts] = ...
            mess_lrradi_get_shifts_hamOpti_generalized( ...
            eqn, opts, oper, W, Z, Y);

    otherwise
        error('MESS:control_data',['unknown shift parameter method: ', ...
            opts.shifts.method]);
end

end
