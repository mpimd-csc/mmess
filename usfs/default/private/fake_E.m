function eqn = fake_E(eqn)

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isfield(eqn,'Ecount'))
    eqn.Ecount = 1;
    eqn.E_ = speye(size(eqn.A_));
else
    eqn.Ecount = eqn.Ecount + 1;
end
