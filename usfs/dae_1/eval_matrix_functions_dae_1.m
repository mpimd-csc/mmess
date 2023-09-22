function [eqn, opts, oper] = eval_matrix_functions_dae_1(eqn, opts, oper, t, sign_dt_E)
%% function eval_matrix_functions_dae_1 updates the matrices in eqn

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if eqn.LTV
    if nargin < 5
        sign_dt_E = 1;
    end

    %%
    if eqn.haveE
        eqn.E_ = eqn.E_time(t);
        eqn.A_ = eqn.A_time(t) + sign_dt_E * eqn.dt_E_time(t);
    else
        eqn.A_ = eqn.A_time(t);
    end
    eqn.B = eqn.B_time(t);
    eqn.C = eqn.C_time(t);
    %% Compute reduced B and C
    n_ode = eqn.manifold_dim;
    if size(eqn.B, 1) > n_ode
        eqn.B = eqn.B(1:n_ode, :) - eqn.A_(1:n_ode, n_ode + 1:end) * ...
                (eqn.A_(n_ode + 1:end, n_ode + 1:end) \ eqn.B(n_ode + 1:end, :));
    end
    if size(eqn.C, 2) > n_ode
        eqn.C = eqn.C(:, 1:n_ode) - ...
                (eqn.C(:, n_ode + 1:end) / ...
                 eqn.A_(n_ode + 1:end, n_ode + 1:end)) * ...
                eqn.A_(n_ode + 1:end, 1:n_ode);
    end
end
end
