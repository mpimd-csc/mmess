function [eqn, opts, oper] = eval_matrix_functions_default(eqn, opts, oper, t, sign_dt_E)
%% function eval_matrix_functions_default updates the matrices in eqn

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
end
end
