function eqn = LU_A(eqn)
% LU_A computes and caches LU decomposition of eqn.A_

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if isfield(eqn, 'AL') && isfield(eqn, 'AU') && ...
        isfield(eqn, 'Ap') && isfield(eqn, 'Aq') && ...
        isfield(eqn, 'AR')
    if isfield(eqn, 'Acount')
        eqn.Acount = eqn.Acount + 1;
    else
        eqn.Acount = 2;
    end
else
    [eqn.AL, eqn.AU, eqn.Ap, eqn.Aq, eqn.AR] = lu(eqn.A_, 'vector');
    eqn.iAq(eqn.Aq) = 1:size(eqn.A_, 1);
    eqn.Acount      = 1;
end
