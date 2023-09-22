function eqn = LU_E_clean(eqn, opts)
% LU_E_CLEAN removes cached LU decomposition of eqn.E_

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if isfield(eqn, 'haveE') && eqn.haveE
    mess_assert(opts, isfield(eqn, 'Ecount'), ...
                'error_arguments', ...
                'field eqn.Ecount is not defined.');

    if eqn.Ecount > 1
        eqn.Ecount = eqn.Ecount - 1;
    else
        eqn = rmfield(eqn, 'Ecount');
        eqn = rmfield(eqn, 'EL');
        eqn = rmfield(eqn, 'EU');
        eqn = rmfield(eqn, 'Ep');
        eqn = rmfield(eqn, 'Eq');
        eqn = rmfield(eqn, 'iEq');
        eqn = rmfield(eqn, 'ER');
    end
end
