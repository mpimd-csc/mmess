function [L, D] = expG(eqn, opts, oper, h, L_zero, D_zero, t_zero)
% Solve the nonlinear problem arising from the split DRE

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

L = L_zero;
Id = eye(size(D_zero));

if not(eqn.LTV) % Autonomous case

    if eqn.type == 'T'
        D = (Id + h * D_zero * L_zero' * ...
             (eqn.B * (eqn.Rinv * (eqn.B' * L_zero)))) \ D_zero;
    elseif eqn.type == 'N'
        D = (Id + h * D_zero * L_zero' * ...
             (eqn.C' * (eqn.Rinv * (eqn.C * L_zero)))) \ D_zero;
    end
else
    % Time-varying case. Only difference is that instead of h*B*Rinv*B' we
    % have int_{t}^{t+h} {B(s)RinvB(s)' ds} we approximate this as
    % LB*DB*LB' by quadrature. Note that the factor h is included in the
    % approximation. For the choice of order 29, see Issue #311.
    [xk, wk] = exact_quadrature_parameters(h, 29);
    xk = xk + t_zero; % Shift from [0, h] t0 [t0, t0+h]

    % Need to evaluate eqn.B or eqn.C' at all xk(k)
    LB = cell(1, length(xk));
    for k = 1:length(xk)
        [eqn, opts, oper] = ...
             opts.splitting.eval_matrix_functions(eqn, opts, oper, xk(k));
        if eqn.type == 'T'
            LB{k} = eqn.B;
        elseif eqn.type == 'N'
            LB{k} = eqn.C';
        end
    end
    LB = cell2mat(LB);

    if eqn.type == 'T'
        DB = kron(diag(wk), eqn.Rinv * speye(size(LB, 2) / length(xk)));
    elseif eqn.type == 'N'
        DB = kron(diag(wk), speye(size(LB, 2) / length(xk)));
    end

    [LB, DB] = mess_column_compression(LB, 'N', DB, ...
                                       opts.splitting.trunc_tol, ...
                                       opts.splitting.trunc_info);

    D = (Id + D_zero * L_zero' * (LB * (DB * (LB' * L_zero)))) \ D_zero;
end
D = mess_symmetrize(D);
end
