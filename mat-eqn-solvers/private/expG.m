function [L, D] = expG(eqn, opts, oper, h, L0, D0, t0)
% Solve the nonlinear problem arising from the split DRE

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

L = L0;
Id = eye(size(D0));

if not(eqn.LTV) % Autonomous case

    if eqn.type == 'T'
        D = (Id + h*D0*L0'*(eqn.B*(eqn.Rinv*(eqn.B'*L0))))\D0;
    elseif eqn.type == 'N'
        D = (Id + h*D0*L0'*(eqn.C'*(eqn.Rinv*(eqn.C*L0))))\D0;
    end
else
    % Time-varying case. Only difference is that instead of h*B*Rinv*B' we
    % have int_{t}^{t+h} {B(s)RinvB(s)' ds} we approximate this as
    % LB*DB*LB' by quadrature. Note that the factor h is included in the
    % approximation. For the choice of order 29, see Issue #311.
    [xk, wk] = exact_quadrature_parameters(h, 29);
    xk = xk + t0; % Shift from [0, h] t0 [t0, t0+h]

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
        DB = kron(diag(wk), eqn.Rinv*speye(size(LB,2) / length(xk)));
    elseif eqn.type == 'N'
        DB = kron(diag(wk), speye(size(LB,2) / length(xk)));
    end

    [LB, DB] = mess_column_compression(LB, 'N', DB, ...
        opts.splitting.trunc_tol, opts.splitting.trunc_info);

    D = (Id + D0*L0'*(LB*(DB*(LB'*L0))))\D0;
end

D = 0.5 * (D + D'); % prevent symmetry loss due to unsymmetric formulas above
end
