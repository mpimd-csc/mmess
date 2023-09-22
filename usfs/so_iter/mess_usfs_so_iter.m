% The second order system
%
%        M*x"(t) + E*x'(t) + K*x(t) = B2*u(t),
%                              y(t) = Cp*x(t) + Cv*x'(t),
%
%  is implicitly transformed to the first order system
%
%                          E_f*x_f' = A_f*x_f + B_f*u,
%                              y(t) = C_f*x_f,
%
%  where
%
%           |  K   αM |
%     E_f = | αM    M | ,
%
%           |-αK   K-αE |
%     A_f = | -K  -E+αM |,
%
%           | αB2|
%     B_f = |  B2|,
%
%     C_f = |Cp Cv|
%
%           | x |
%     x_f = | x'|,
%
% and 0 < α < λ_min( E / (M+.25E*K\E) ), following [1] (Thm 6.1).
%
% The Matrix M, E, K are assumed to be symmetric positive definite and square.
% The fieldnames have to end with _ to indicate that the data are
% inputdata for the Algorithm.
%
% eqn.M_ = M
% eqn.K_ = K
% eqn.E_ = E
% eqn.B  = B_f
% eqn.C  = C_f
%
% In contrast to the other second order usfs, here all solve operations are
% performed using iterative solvers rather than "\".
% The `sol_*` calls will look into `opts.usfs.so_iter` for
% further setting. `opts.usfs.so_iter` is a structure with
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
% Moreover, `opts.usfs.so_iter` should contain
%
%    alpha       α in the implicitly transformed system above
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
% References:
%
% [1] H. K. F. Panzer, Model order reduction by Krylov subspace
%     methods with global error bounds and automatic choice of
%     parameters, Dissertation, Technische Universität München,
%     Munich, Germany (2014).
%     https://mediatum.ub.tum.de/doc/1207822/1207822.pdf
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
