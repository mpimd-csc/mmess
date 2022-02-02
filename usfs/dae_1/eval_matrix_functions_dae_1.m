function [ eqn, opts, oper ] = eval_matrix_functions_dae_1( eqn, opts, oper, t )
%% function eval_matrix_functions_dae_1 updates the matrices in eqn

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if eqn.LTV

    %%
    if eqn.haveE
        eqn.E_ = eqn.E_time(t);
        eqn.A_ = eqn.A_time(t) + eqn.dt_E_time(t);
    else
        eqn.A_ = eqn.A_time(t);
    end
    eqn.B = eqn.B_time(t);
    eqn.C = eqn.C_time(t);
    %% Compute reduced B and C
    st = eqn.st;
    if size(eqn.B, 1) > st
        eqn.B = eqn.B(1 : st, :) - eqn.A_(1 : st, st + 1 : end) ...
            * (eqn.A_(st + 1 : end, st + 1 : end) \ eqn.B(st + 1 : end, :));
    end
    if size(eqn.C, 2) > st
        eqn.C = eqn.C( : , 1 : st) - (eqn.C( : , st + 1 : end) ...
            / eqn.A_(st +1 : end, st + 1 : end)) * eqn.A_(st+1 : end, 1 : st);
    end
end
end