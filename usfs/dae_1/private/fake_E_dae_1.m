function eqn = fake_E_dae_1(eqn)
% FAKE_E_DAE_1 adds a dummy eqn.E_, for later reference in ApE
% routines, when "haveE" is 'false'.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isfield(eqn, 'Ecount'))
    eqn.Ecount = 1;
    % E = [ I 0 ]
    %     [ 0 0 ]
    eqn.E_ = sparse(1:eqn.manifold_dim, ...
                    1:eqn.manifold_dim, ...
                    ones(eqn.manifold_dim, 1), ...
                    n, n, eqn.manifold_dim);
else
    eqn.Ecount = eqn.Ecount + 1;
end
