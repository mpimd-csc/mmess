% The first order system
%
%      E * z'(t) = A * z(t) + B * u(t)
%           y(t) = C * z(t)
%
% is encoded in the eqn structure
%
% The fieldnames for A and E have to end with _  to indicate that the data
% are inputdata for the algorithm.
%
% eqn.A_ = A
% eqn.E_ = E
% eqn.B  = B
% eqn.C  = C
%
% Note that E_ and A_ are expected to be sparse and of size n x n,
% while B and C may be dense (and will be converted to dense by
% some routines anyway) and should have far less columns (for B)
% and rows (for C) than n.
%
% In contrast to the "default" usfs, here all solve operations are
% performed using iterative solvers rather than "\".
% The `sol_*` calls will look into `opts.usfs.default_iter` for
% further setting. `opts.usfs.default_iter` is a structure with
% members:
%
%    method_A     the iterative solver function for solving with A
%
%    method_ApE   the iterative solver function for the shifted solves
%
%    method_E     the iterative solver function for solving with E
%
% all of these can currently be any suitable iterative solver from
% MATLAB or GNU Octave. we are regularly testing:
%
% matlab_solvers = {'bicg', 'bicgstab', 'bicgstabl', ...
%    'cgs', 'gmres', 'lsqr', 'minres', 'pcg', 'qmr', 'symmlq', 'tfqmr'};
% octave_solvers = {'bicg', 'bicgstab', ...
%    'cgs', 'gmres', 'pcg', 'qmr', 'tfqmr', 'pcr'};
%
% for `method_A` and `method_E`, while fixing `method_ApE = 'gmres'`
%
% Moreover, `opts.usfs.default_iter` should contain
%
%    res_tol     residual tolerance passed to all iterative solvers
%
%    max_iter    maximum iteration number passed to all iterative solvers
%
%    restIter    restart length passed to GMRES
%
%    PA_L, PA_R  preconditioner matrices for A passed to the iterative
%                solvers (as M1 and M2 in the default MATLAB
%                iterative solver interface)
%    PE_L, PE_R  preconditioner matrices for E passed to the iterative
%                solvers (as M1 and M2 in the default MATLAB
%                iterative solver interface)
%
%  `oper.init` will check, and, in case they are absent, initialize them
%  with defaults.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
