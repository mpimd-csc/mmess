function [rw, Hp, Hm, Vp, Vm, eqn, opts, oper] = ...
    get_ritz_vals_dae_3_so(eqn, opts, oper, U, W, p_old)
%  This function is an exact copy of the Penzl heuristic part in mess_para.
%  the only difference is that B or C and K are filled up by trailing zero
%  blocks to allow for the computation of the Ritz values with respect to
%  the full size block matrices instead of the restriction to the
%  (1,1)-block for which the ADI is formulated.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% Input data not completely checked!
for mat = 'MEKG'
    if not(isfield(eqn, sprintf('%c_', mat))) || ...
       not(eval(sprintf('isnumeric(eqn.%c_)', mat)))
        mess_err(opts, 'error_arguments', 'field eqn.%c_ is not defined', mat);
    end
end

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');
if not(result)
    mess_err(opts, 'control_data', ...
             'system data is not completely defined or corrupted');
end
% returns order of A or states of A, A is supposed to be square
nv = size(eqn.M_, 1);
np = size(eqn.G_, 1);

%%
% here we add the trailing zero blocks. Note that we are passing the
% eqn structure back as an output, so to ensure this change is not visible in
% anything above this routine and will only be passed on to the function
% handles used in here, we need to truncate again later.
if isfield(eqn, 'U') && not(isempty(eqn.U))
    eqn.U = [eqn.U; sparse(np, size(eqn.U, 2))];
end
if isfield(eqn, 'V') && not(isempty(eqn.V))
    eqn.V = [eqn.V; sparse(np, size(eqn.V, 2))];
end

if isfield(opts.shifts, 'method') && ...
        strcmp(opts.shifts.method, 'projection')
    U = [U; zeros(2 * nv + np - size(U, 1), size(U, 2))];
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
    else
        W = [W; zeros(2 * nv + np - size(W, 1), size(W, 2))];
    end
    rw = mess_projection_shifts(eqn, opts, oper, U, W, p_old);
else
    if not(isfield(opts.shifts, 'b0')) || isempty(opts.shifts.b0)
        opts.shifts.b0 = ones(n, 1);
    else
        if not(length(opts.shifts.b0) == 2 * nv + np)
            mess_warn(opts, 'b0', ...
                      'b0 has the wrong length. Switching to default.');
            opts.shifts.b0 = ones(2 * nv + np, 1);
        end
    end
    [rw, Hp, Hm, Vp, Vm] = mess_get_ritz_vals(eqn, opts, oper);
end
if isfield(opts.shifts, 'truncate') && isnumeric(opts.shifts.truncate)
    rw = rw(abs(rw) < opts.shifts.truncate);
    rw = rw(abs(rw) > 1 / opts.shifts.truncate);
end
%%
% Let's truncate U and V back
if isfield(eqn, 'U') && not(isempty(eqn.U))
    eqn.U = eqn.U(1:2 * nv, :);
end
if isfield(eqn, 'V') && not(isempty(eqn.V))
    eqn.V = eqn.V(1:2 * nv, :);
end

end
