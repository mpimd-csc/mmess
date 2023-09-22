function [eqn, opts, oper] = sol_ApE_post_so_iter(eqn, opts, oper)
% function [eqn, opts, oper] = sol_ApE_post_so_iter(eqn, opts, oper)
% It is necessary to remove the identity added in _pre_ again

% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)

if not(eqn.haveE)
    eqn = fake_E_clean_default(eqn);
end
