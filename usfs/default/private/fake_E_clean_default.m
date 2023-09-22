function eqn = fake_E_clean_default(eqn)

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if eqn.Ecount > 1
    eqn.Ecount = eqn.Ecount - 1;
else
    eqn = rmfield(eqn, {'E_', 'Ecount'});
end
