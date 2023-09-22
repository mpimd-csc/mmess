function eqn = LU_A_clean(eqn, opts)
% LU_A_CLEAN removes cached LU decomposition of eqn.A_

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
mess_assert(opts, isfield(eqn, 'Acount'), ...
            'error_arguments', ...
            'field eqn.Acount is not defined.');

if eqn.Acount > 1
    eqn.Acount = eqn.Acount - 1;
else
    eqn = rmfield(eqn, 'Acount');
    eqn = rmfield(eqn, 'AL');
    eqn = rmfield(eqn, 'AU');
    eqn = rmfield(eqn, 'Ap');
    eqn = rmfield(eqn, 'Aq');
    eqn = rmfield(eqn, 'iAq');
    eqn = rmfield(eqn, 'AR');
end
