function [eqn, opts, oper] = mul_ApE_post_dae_1(eqn, opts, oper)
% we need to remove the identity added in _pre_ again.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(eqn.haveE)
    eqn = fake_E_clean(eqn);
end