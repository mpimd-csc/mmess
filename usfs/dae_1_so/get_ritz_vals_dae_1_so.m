function [rw, Hp, Hm, Vp, Vm, eqn, opts, oper] = ...
    get_ritz_vals_dae_1_so(eqn, opts, oper, U, W, p_old)
% [rw, Hp, Hm, Vp, Vm, eqn, opts, oper] = ...
%    get_ritz_vals_dae_1_so(eqn, opts, oper, U, W, p_old)
%
% Wrapper for the special system structure around mess_get_ritz_vals.
% Additionally due to the second order structure, the real value
% opts.shifts.truncate can be set to remove any computed values that are
% smaller than opts.shifts.truncate, or larger than 1/opts.shifts.truncate.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% Input data not completely checked!
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');
if not(result)
    mess_err(opts, 'control_data', ...
             'system data is not completely defined or corrupted');
end
n = oper.size(eqn, opts);

%%

if isfield(opts.shifts, 'method') && ...
        strcmp(opts.shifts.method, 'projection')
    if isempty(W)
        % first shifts are computed with U = eqn.W and W = A * eqn.W
        W = oper.mul_A(eqn, opts, eqn.type, U, 'N');
        if isfield(eqn, 'haveUV') && eqn.haveUV
            switch eqn.type
                case 'N'
                    W = W + eqn.U * (eqn.V' * U);
                case 'T'
                    W = W + eqn.V * (eqn.U' * U);
            end
        end
    end
    rw = mess_projection_shifts(eqn, opts, oper, U, W, p_old);
else
    if not(isfield(opts.shifts, 'b0')) || isempty(opts.shifts.b0)
        opts.shifts.b0 = ones(n, 1);
    else
        if not(length(opts.shifts.b0) == n)
            mess_warn(opts, 'b0', ...
                      'b0 has the wrong length. Switching to default.');
            opts.shifts.b0 = ones(n, 1);
        end
    end
    [rw, Hp, Hm, Vp, Vm] = mess_get_ritz_vals(eqn, opts, oper);
end
%
% remove Ritz values that are too large or too small if desired
if isfield(opts.shifts, 'truncate') && isnumeric(opts.shifts.truncate)
    rw = rw(abs(rw) < opts.shifts.truncate);
    rw = rw(abs(rw) > 1 / opts.shifts.truncate);
end
