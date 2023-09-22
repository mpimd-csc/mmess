function [opts, l] = mess_get_projection_shifts(eqn, opts, oper, Z, W, D)
% Intended for use inside the ADI iteration,
% mess_get_projection_shifts  computes new projection shifts and
% updates shift vectors if the shift computation method is
% 'projection', otherwise it simply does nothing. The function is
% called in mess_lradi whenever the end of the current shift vector
% is reached.
%
% Inputs:
%  eqn, opts, oper  the usual structures containing the equation
%                   data, the process control parameters and the
%                   operation function handles.
%
%  Z, W, D          The matrices for the current residual (W) and
%                   the tall and skinny factor in the solution
%                   approximation (Z) and in case of LDL_T D the
%                   factorization kernel
%
% Output:
%  opts             altered opts structure, with the new shift vector in
%                   opts.shifts.p
%  l                the number of admissible shifts computed,
%                   i.e. the length of opts.shifts.p
%
%  mess_get_projection shifts evaluates opts.shifts.num_desired to
%  determine how many shifts are to be computed and how many block
%  columns from Z need to be passed to mess_projection_shifts to
%  compute them.
%
%  By default a shift is admissible if it resides in the left half
%  of the complex plane.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check input
ncols_W = size(eqn.W, 2);

if not(isfield(opts, 'shifts')) || not(isstruct(opts.shifts))
    mess_warn(opts, 'control_data', ...
              ['shift parameter control structure missing. ', ...
               'Switching to default num_desired = 25.']);
    opts.shifts.num_desired = 25;
else
    if not(isfield(opts.shifts, 'num_desired')) || ...
            not(isnumeric(opts.shifts.num_desired))
        mess_warn(opts, 'control_data', ...
                  ['Missing or Corrupted opts.shifts.num_desired field. ', ...
                   'Switching to default: 25']);
        opts.shifts.num_desired = 25;
    end
end

if not(isfield(opts.shifts, 'banned')) || ...
        not(isnumeric(opts.shifts.banned))
    opts.shifts.banned = [];
elseif not(isfield(opts.shifts, 'banned_tol')) || ...
        not(isnumeric(opts.shifts.banned_tol)) || ...
        not(isscalar(opts.shifts.banned_tol))
    opts.shifts.banned_tol = 1e-4;
end

if not(isfield(opts.shifts, 'recursion_level')) || ...
        not(isnumeric(opts.shifts.recursion_level)) || ...
        not(isscalar(opts.shifts.recursion_level))
    opts.shifts.recursion_level = 0;
end

if isfield(opts.shifts, 'info') && opts.shifts.info > 0
    info = 1;
else
    info = 0;
end
%% Compute new shifts
if isfield(opts.shifts, 'method') && ...
        strcmp(opts.shifts.method, 'projection')

    if info
        mess_fprintf(opts, 'updating shifts\n');
    end

    if not(isfield(opts.shifts, 'used_shifts')) || ...
            isempty(opts.shifts.used_shifts)
        opts.shifts.used_shifts = opts.shifts.p;
    else
        opts.shifts.used_shifts = ...
            [opts.shifts.used_shifts; opts.shifts.p];
        first_dropped = length(opts.shifts.used_shifts) - ...
                        opts.shifts.num_desired + 1;
        last_kept = first_dropped - 1;
        if (size(opts.shifts.used_shifts, 1) > opts.shifts.num_desired) && ...
                imag(opts.shifts.used_shifts(first_dropped)) && ...
                (abs(opts.shifts.used_shifts(first_dropped) - ...
                     conj(opts.shifts.used_shifts(last_kept))) < eps)
            % don't cut between pair of complex shifts
            opts.shifts.used_shifts = ...
                opts.shifts.used_shifts(end - opts.shifts.num_desired:end);
        elseif size(opts.shifts.used_shifts, 1) > opts.shifts.num_desired
            opts.shifts.used_shifts = ...
                opts.shifts.used_shifts(end - opts.shifts.num_desired + 1:end);
        end
    end

    %% Compute new shifts
    if isfield(oper, 'get_ritz_vals')
        if opts.LDL_T
            % scale columns of Z (L) as in original non LDL^T formulation
            len = size(opts.shifts.used_shifts, 1) * ncols_W - 1;
            inds_D = size(D, 1) - size(opts.shifts.used_shifts, 1) + ...
                     1:size(D, 1);
            p = oper.get_ritz_vals(eqn, opts, oper, ...
                                   Z(:, end - len:end) * ...
                                   kron(sqrt(D(inds_D, inds_D)), ...
                                        eye(ncols_W)), ...
                                   W, opts.shifts.used_shifts);
        else
            len = (size(opts.shifts.used_shifts, 1) * ncols_W) - 1;
            p = oper.get_ritz_vals(eqn, opts, oper, ...
                                   Z(:, end - len:end), ...
                                   W, opts.shifts.used_shifts);
        end
    else
        if opts.LDL_T
            % scale columns of Z (L) as in original non LDL^T formulation
            len = size(opts.shifts.used_shifts, 1) * ncols_W - 1;
            inds_D = size(D, 1) - size(opts.shifts.used_shifts, 1) + ...
                     1:size(D, 1);
            p = mess_projection_shifts(eqn, opts, oper, ...
                                       Z(:, end - len:end) * ...
                                       kron(sqrt(D(inds_D, inds_D)), ...
                                            eye(ncols_W)), ...
                                       W, opts.shifts.used_shifts);
        else
            len = size(opts.shifts.used_shifts, 1) * ncols_W - 1;
            p = mess_projection_shifts(eqn, opts, oper, ...
                                       Z(:, end - len:end), ...
                                       W, opts.shifts.used_shifts);
        end
    end

    %% check computed shifts
    % check for banned shifts
    for j = 1:length(opts.shifts.banned)
        critical_shifts = abs(p - opts.shifts.banned(j)) < ...
            opts.shifts.banned_tol * max(abs(p));
        p(critical_shifts) = p(critical_shifts) - ...
            opts.shifts.banned_tol * max(abs(p)) * 2;
    end
    if isempty(p) % if all shifts banned try again with double amount
        if opts.shifts.recursion_level < 2
            mess_warn(opts, 'projection_shifts', ...
                      'All computed shifts have been banned. Retrying');
            num_desired = opts.shifts.num_desired;
            opts.shifts.num_desired = num_desired * 2;
            opts.shifts.used_shifts = ...
                opts.shifts.used_shifts(1:end - size(opts.shifts.p, 1));
            opts.shifts.recursion_level = opts.shifts.recursion_level + 1;
            if opts.LDL_T
                [opts] = mess_get_projection_shifts(eqn, opts, oper, Z, W, D);
            else
                [opts] = mess_get_projection_shifts(eqn, opts, oper, Z, W);
            end
            opts.shifts.num_desired = num_desired;
            p = opts.shifts.p;
            opts.shifts.recursion_level = opts.shifts.recursion_level - 1;
        end
    end
    if not(isempty(p))
        opts.shifts.p = p;

        if info
            mess_fprintf(opts, 'p:\n');
            for ip = 1:length(p)
                mess_fprintf(opts, '%e\n', p(ip));
            end
        end
    else
        % could not compute new shifts, reuse previous ones
        mess_warn(opts, 'projection_shifts', ...
                  'projection update returned empty set. ', ...
                  'Reusing previous set!');
    end
end
l = length(opts.shifts.p);
