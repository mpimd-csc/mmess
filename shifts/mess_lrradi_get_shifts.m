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
% opts.shifts.method    possible  values:
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
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isfield(opts.shifts, 'banned')) || ...
        not(isnumeric(opts.shifts.banned))
    opts.shifts.banned = [];
elseif not(isfield(opts.shifts, 'banned_tol')) || ...
        not(isnumeric(opts.shifts.banned_tol)) || ...
        not(isscalar(opts.shifts.banned_tol))
    opts.shifts.banned_tol = 1e-4;
end

%% Check input
switch opts.shifts.method
    case 'precomputed'
        % Just return the same array of shifts as the user provided.
        nShifts = length(opts.shifts.p);

    case 'heur'
        % Use MESS routines for heuristic penzl shifts.
        [p, ~, eqn, opts, oper] = mess_para(eqn, opts, oper);
        nShifts                 = length(p);
        opts.shifts.p           = p;

    case 'projection'
        % Use MESS routines for projection shifts.
        if isempty(Z)
            Z = W;
        end

        nShifts = 0;
        i = 1;
        while nShifts == 0
            nZ = size(Z, 2);
            if not(opts.radi.compute_sol_facpart)
                % Relevant parts of LRF.
                maxcolZ = opts.shifts.history;
                ind     = max(nZ - maxcolZ, 0) + 1:nZ;
            else
                ind = 1:nZ;
            end

            if eqn.type == 'T'
                p = mess_projection_shifts(eqn, opts, oper, Z(:, ind), ...
                                           oper.mul_A(eqn, opts, 'T', ...
                                                      Z(:, ind), 'N') + ...
                                           eqn.V * (eqn.U' * Z(:, ind)), []);
            else
                p = mess_projection_shifts(eqn, opts, oper, Z(:, ind), ...
                                           oper.mul_A(eqn, opts, 'N', ...
                                                      Z(:, ind), 'N') + ...
                                           eqn.U * (eqn.V' * Z(:, ind)), []);
            end

            for j = 1:length(opts.shifts.banned)
                critical_shifts = abs(p - opts.shifts.banned(j)) < ...
                    opts.shifts.banned_tol * max(abs(p));
                p(critical_shifts) = p(critical_shifts) - ...
                    opts.shifts.banned_tol * max(abs(p)) * 2;
            end

            nShifts       = length(p);
            opts.shifts.p = p;
            if isempty(p)
                if i < 5
                    mess_warn(opts, 'mess_para', ...
                              ['Could not compute initial projection ', ...
                               'shifts. Going to retry with random ' ...
                               'right hand side.']);
                    Z = rand(size(Z));
                else
                    mess_err(opts, 'mess_para', ...
                             'Could not compute initial projection shifts.');
                end
            end
            i = i + 1;
        end
    case 'gen-ham-opti'
        % Residual hamiltonian shifts... recommended+default!
        [eqn, opts, oper, nShifts] = ...
            mess_lrradi_get_shifts_hamOpti_generalized(eqn, opts, oper, ...
                                                       W, Z, Y);

    otherwise
        mess_err(opts, 'control_data', ...
                 ['unknown shift parameter method: ', ...
                  opts.shifts.method]);
end

end
