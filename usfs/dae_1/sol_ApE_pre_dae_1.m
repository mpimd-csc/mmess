function [eqn, opts, oper] = sol_ApE_pre_dae_1(eqn, opts, oper)
% to simplify matters in sol_ApE we add a field eqn.E_ holding the
% identity matrix when we do not have an E matrix already.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(eqn.haveE)
    eqn = fake_E_dae_1(eqn);
end
