%
%  ######################################################################
%  #                                                                    #
%  #          M-M.E.S.S. - Matrix Equations, Sparse Solvers             #
%  #                                                                    #
%  ######################################################################
%
%  version 3.0
%
%  The M-M.E.S.S. toolbox is intended for the solution of symmetric
%  linear and quadratic, differential and algebraic matrix
%  equations with real, large, and sparse (or sparse + low-rank)
%  coefficient matrices in the linear terms, low-rank quadratic and
%  constant terms, and low-rank solutions.
%
%  Using the user supplied functions (usfs) framework,
%  M-M.E.S.S. enables the solution of equations for several
%  system structures. While implicitly M-M.E.S.S. always treats a
%  system
%      .
%    E x(t) = A x(t) + B u(t),                                 (1)
%      y(t) = C x(t),
%
%  with E invertible, the usfs support certain structured
%  differential algebraic equation (DAE) systems, first-order forms
%  of second-order differential equations and combinations of
%  both. (See section USFS below)
%
%  Features:
%   * large-scale algebraic Lyapunov equations
%     (mess_lyap, mess_lradi, mess_KSM)
%   * large-scale algebraic Lyapunov-plus-positive equations
%     (mess_lyapunov_bilinear)
%   * large-scale algebraic Riccati equations
%     (mess_care, mess_lrnm, mess_lrri, mess_lrradi, mess_KSM)
%   * large-scale differential Riccati equations
%     (mess_bdf_dre, mess_rosenbrock_dre, mess_splitting_dre)
%   * model order reduction
%     (mess_balanced_truncation, mess_square_root_method,
%      mess_tangential_IRKA)
%
%  Moreover a first solver for certain, so called sparse-dense Sylvester
%  equations has been added in mess_sylvester_sparse_dense and
%  Lyapunov-plus-positive equations related to bilinear control systems
%      .
%    E x(t) = A x(t) + sum_k N_k u_k x + B u(t),               (2)
%      y(t) = C x(t),
%
%  are supported.
%
%  ######################################################################
%  #                                                                    #
%  #  supported matrix equations                                        #
%  #                                                                    #
%  ######################################################################
%
%  1.) algebraic Lyapunov equations
%
%    A X E' + E X A' + B B' = 0
%    A' X E + E' X A + C' C = 0
%
%  For the case where E, A are real sparse matrices we provide a lyap-like
%  function mess_lyap. In general these equations are solved via the
%  low-rank (X = Z Z', or X = L D L') alternating directions implicit
%  iteration (LR-ADI) using the function mess_lradi. Note that in contrast
%  to the lyap case, our functions work with factorized constant terms and
%  solutions, i.e. the mess_lyap call would be:
%
%    Z = mess_lyap(A, B, [], [], E)
%
%  for the first equation above. Alternatively, these equations can be
%  solved by mess_KSM for classic first-order systems.
%
%  2.) continuous time algebraic Riccati equations
%
%    A X E' + E X A' - E X B B' X E' + C' C = 0
%    A' X E + E' X A - E' X C' C X E + B B' = 0
%
%  As in the Lyapunov case, the solution can be approximated in Z Z' and
%  L D L' format.
%  For these control and filter equations, with respect to standard systems
%  of the form (1), we provide several solvers. Similar to mess_lyap above
%  we also provide a quick-start interface mess_care using a similar syntax
%  as care, i.e.
%
%   [Z, K] = mess_care(A, B, C, [], E)
%
%  computes the solution factor Z such that we approximately have X = Z Z'
%  and analogous to care K is the feedback matrix stabilizing (1).
%
%  For experts and for use with the non-standard usfs such that system (1)
%  is only given implicitly, we implemented a low-rank Newton ADI method
%  that supports inexact Newton iterations and line search in
%  mess_lrnm, as well as the RADI iteration in mess_lrradi. While mess_lrnm
%  can compute the solution in the above formats and supports computing the
%  feedback K without ever forming Z, mess_lrradi computes the solution as
%  L inv(D) L'. Alternatively, these equations can be solved by mess_KSM
%  for classic first-order systems.
%
%  3.) Riccati equations with indefinite quadratic terms
%
%    A X E' + E X A' + E X (C1' C1 - C2' C2) X E' + B1 B1' = 0
%    A' X E + E' X A + E' X (B1' B1 - B2' B2) X E + C1' C1 = 0
%
%  For these control and filter equations (arising in H-infinity control)
%  we provide a low-rank Riccati iteration solver (X = Z Z').
%  Especially, this covers the standard continuous-time algebraic Riccati
%  equation (see 2.)) and the positive algebraic Riccati equations
%
%    A X E' + E X A' + E X C1' C1 X E' + B1 B1' = 0,
%    A' X E + E' X A + E' X B1' B1 X E + C1' C1 = 0.
%
%  4.) differential Riccati equations
%                                                .
%    A X E' + E X A' - E X B B' X E' + C' C = -E X E'
%                                                 .
%    A' X E + E' X A - E' X C' C X E + B B' = -E' X E
%
%  We support backward differentiation formulae as 1 to 4 step
%  methods or Rosenbrock 1-stage and 2-stage schemes in the
%  functions mess_bdf_dre and mess_rosenbrock_dre. In principle both
%  solvers can also solve differential Lyapunov equations by
%  setting the data for the quadratic term to zero. Specialized
%  more efficient solvers exploiting this have not yet been
%  implemented. Note that these solvers necessarily use the form
%
%    X = L D L'.
%
%  The differential Riccati solver support both the autonomous case (i.e.
%  (1) is linear time-invariant) as well as the non-autonomous case (i.e.
%  (1) is a linear time-varying system).
%
%  5.) sparse-dense Sylvester equations
%
%   A * X * F + E * X * H = -M
%
%  where A, E are the large and sparse matrices from (1) and F, H are
%  usually small to moderate size dense matrices, such that M and X are
%  tall rectangular.
%
%  6.) algebraic Lyapunov-plus-positive equations
%
%   A*P*E' + E*P*A' + Sum_N_k*P*N_k' + B*B' = 0       (N)
%   A'*Q*E + E'*Q*A + Sum_N_k'*Q*N_k + C'*C = 0       (T)
%
%  where in addition to the usual coefficients the sparse and square matrices
%  N_k (number of columns in B many in (N) or rows in C many in (T)) enter
%  the problem. Currently this equation is only supported by the default
%  usfs.
%
%  ######################################################################
%  #                                                                    #
%  #  supported model reduction methods                                 #
%  #                                                                    #
%  ######################################################################
%
%  Currently Balanced Truncation (BT) and (tangential) IRKA are supported
%  for systems in explicit form (1) in the routines
%  mess_balanced_truncation and mess_tangential_irka. For (1) given by its
%  coefficient matrices (and E invertible) a reduced order model of the
%  same form can be computed by
%
%  [Er, Ar, Br, Cr, outinfo] =
%     mess_balanced_truncation(E, A, B, C, [], opts)
%
%  or
%
%  [Er, Ar, Br, Cr, S, b, c, V, W] = mess_tangential_irka(E, A, B, C, opts)
%
%  where opts.bt.max_order, opts.bt. trunc_tol set the maximum desired
%  reduced order and the truncation tolerance for the BT error bound, while
%  opts.irka.r, opts.irka.maxiter, opts.irka.shift_tol, opts.irka.h2_tol
%  stand for the reduced order, them maximum allowed IRKA steps, the tolerance
%  for the relative change of shifts stopping criterion and the relative
%  change of the H2 system norm for subsequent admissible iterates in IRKA.
%
%  Several demonstration functions in the DEMOS sub-folders show how BT can
%  be performed for other system structures implemented by the USFS below.
%
%  Furthermore, helper tasks are implemented in
%  mess_square_root_method, mess_tf_plot. Note that with release 3.0 these
%  functions have been deprecated and will not be maintained or developed
%  any longer. Instead the development is shifted to the MORLAB package
%  that will use M-M.E.S.S. as the solver backend for the required matrix
%  equations for model order reduction of sparse systems starting from
%  version 6.0.
%
%  ######################################################################
%  #                                                                    #
%  #  USFS - user supplied functions and system structures              #
%  #                                                                    #
%  ######################################################################
%
%  All operations with the implicit matrices E, A, and (A + p E) regarding
%  multiplication and solution of linear systems are implemented via user
%  supplied functions used via function handles in the oper structure. The
%  following describes the sets of function handles shipped with M-M.E.S.S.
%  and what types of systems they implement.
%
%  1.) generalized first order systems (1): "default", "default_iter",
%      "state_space_transformed_default"
%    The default set of of function handles is "default". These function handles
%    directly act on the matrices E, A given in (1) using them to
%    explicitly represent the actions. That means multiplications directly
%    use A* and E* (or their transposes if requested) and linear solves are
%    implemented via A\, E\ (A + p*E)\. Again with transposes where
%    requested. For the "default_iter" set all \ calls are replaced by
%    configurable iterative linear solvers. The
%    "state_space_transformed_default" set uses \ but transforms the system
%    to standard statespace form (E=I) for compatibility with the standard
%    Krylov projection method formulations.
%
%  2.) second order: "so_1", "so_2"
%    These sets of function handles work with second order dynamical
%    systems
%         ..        .
%      M  z (t) + E z(t) + K z(t) = B2 u(t)
%                .
%      y(t) = Cv z(t) + Cp z(t)
%                                                              .
%    using their representation in phase space, i.e. x = [ z ; z ], with
%    different representations of the coefficients E and A in (1).
%    See the Readme files in the corresponding folders for details on the
%    phase space representations. Note that these are only used implicitly
%    to formulate the algorithm in therms of (1) , but all operations are
%    performed in terms of the original matrices M, E, K, i.e. in n instead
%    of 2n degrees of freedom.
%
%
%  3.) semi-explicit index-1 DAEs: "dae_1", "dae_1_so"
%    The "dae_1" set of usfs expects
%
%      Ef = [ E1 0 ]    Af = [ A1 A2 ]
%           [ 0  0 ]         [ A3 A4 ]
%
%    and A4 invertible. It then implicitly computes with E = E1 and
%    A = A1 - A2*A4\A3, but never forms these explicitly. The "dae_1_so"
%    set combines this with phase space representation of second order
%    systems.
%
%  4.) Stokes-type index-2 DAEs: "dae_2"
%    The "dae_2" set of USFS is intended for proper index-2 DAE systems of
%    the structure
%
%      Ef = [ E1 0 ]    Af = [ A1 G ]
%           [ 0  0 ]         [ H  0 ]
%
%    and implicitly computes with the matrices restricted to the hidden
%    manifold on which (1) represents an ODE. This representation is never
%    formed explicitly, but some of the equation data, e.g. the factors of
%    the constant terms, are explicitly projected (requiring additional
%    memory).
%
%  5.) index-2 and index-3 constraint mechanical systems: "dae_2_so",
%      "dae_3_so"
%    After phase space representation as in 2.) above, second order systems
%    with constraints purely on the positions (index-3) or velocities
%    (index-2) immediately take the form in 4.), such that again one can
%    implicitly compute on the hidden manifold, without forming the
%    projected coefficients, but only projecting certain data and avoiding
%    the doubling of dimensions with the analogous techniques as in 2.).
%
% For further details see also:
%
%    help mess_usfs
%    help mess_usfs_dae_1
%    help mess_usfs_dae_1_so
%    help mess_usfs_dae_2
%    help mess_usfs_default_iter
%    help mess_usfs_default
%    help mess_usfs_so_1
%    help mess_usfs_so_2
%    help mess_usfs_so_iter
%    help mess_usfs_state_space_transformed_default
%
%  ######################################################################
%  #                                                                    #
%  #  SPLR operators - A in (1) is sparse with a low-rank update        #
%  #                                                                    #
%  ######################################################################
%
%  We consider the form
%
%    A = F + U V'
%
%  for a sparse matrix F and tall and skinny rectangular matrices U, V and
%  exploit the Sherman-Morrison-Woodbury formula
%              -1    -1    -1             -1    -1     -1
%    (F + U V')   = F   - F   U ( I + V' F   U )   V' F
%
%  for computing the action of the inverse of A. Due to the USFS setup this
%  works for all implicit equations (1). Note further that we have
%  developed the codes with U V' implementing some sort of stabilization in
%  mind such that we believe the usual numerical stability doubts about
%  this formula do not apply, here.
%
%  ######################################################################
%  #                                                                    #
%  #  eqn, opts, oper                                                   #
%  #                                                                    #
%  ######################################################################
%
%  All mess routines other than the easy wrapper interfaces mess_lyap and
%  mess_care, as well as the two MOR routines mess_balanced_truncation and
%  mess_tangential_irka depend on the three default structures eqn, opts,
%  oper. Their general purposes are described in the following. The
%  specific content depends on the methods used, and the system structure
%  handled by the chosen USFS. All three fields are input and output
%  parameters to basically all routines and can thus also be used to store
%  additional user data that may be needed by USFS not included with
%  M-M.E.S.S. A prototypical mess function call would thus look like
%
%  [*, eqn, opts, oper] = mess_xxxx(eqn, opts, oper, *)
%
%  eqn
%  ---
%  is a struct that stores the user supplied problem data describing (1).
%  In the simplest case it just contains the matrices E, A, B, C.
%
%  opts
%  ----
%  hosts all information used to control the algorithms. It selects the
%  format of the desired solution, the norms for residual norm evaluations
%  in stopping criteria, tolerances, flags steering variants of the
%  algorithms and alike. It generally consists of a number of substructures
%  per algorithm involved in the program. The members and names of these
%  structures are given in the help texts of the single M-M.E.S.S.
%  routines.
%
%  oper
%  ----
%  stores all information regarding the USFS. It holds references to all
%  user supplied functions and may be used to store additional data
%  required by these functions.
%
%  ######################################################################
%  #                                                                    #
%  #  Logging                                                           #
%  #                                                                    #
%  ######################################################################
%
%  The M-M.E.S.S. routines support a custom logging mechanism. This allows
%  the user to select the output of the routines to be rerouted to files,
%  the command window or both. As an effect, mess no longer uses the
%  standard warnings and errors, but the following custom functions:
%
%  mess_log_initialize
%  ----
%  initializes the logging with user set parameters
%
%  mess_log_finalize
%  ----
%  ends logging and properly closes all related files
%
%  mess_fprintf
%  ----
%  prints the given string to desired output,
%  accepts sprintf like syntax
%
%  mess_err
%  ----
%  issuing errors to the desired output
%
%  mess_warn
%  ----
%  issuing warnings to the desired output
%
%  mess_assert
%  ----
%  assertions sending error messages to the desired output
%
%  mess_log_plot
%  ----
%  saves figures in the set format and includes them into the
%  log document if specified
%
%  mess_log_matrix
%  ----
%  saves the given variable to a .mat file and issues a notification
%  in the logging output
%
%  For more insight on the logging output, refer to DEMOS/logging.m for a
%  demo program
%
%
%  ######################################################################
%  #                                                                    #
%  #  Citation                                                          #
%  #                                                                    #
%  ######################################################################
%
%  If you want to refer to the software in general please cite:
%
%   J. Saak, M. Koehler, P. Benner, M-M.E.S.S. – the matrix equations sparse
%   solvers library,
%   see also: https://www.mpi-magdeburg.mpg.de/projects/mess.
%   DOI: 10.5281/zenodo.632897.
%
%  In case you want to refer to this specific version please check the
%  CITATION.md file in the installation folder.
%
%  For the theoretic backing and basic philosophy we recommend citing
%
%   P. Benner, M. Koehler, J. Saak, Matrix equations, sparse solvers:
%   M-M.E.S.S.-2.0.1 – philosophy, features and application for (parametric)
%   model order reduction,
%   in: P. Benner, T. Breiten, H. Faßbender, M. Hinze, T. Stykel, R. Zimmermann
%   (Eds.), *Model Reduction of Complex Dynamical Systems*, Vol. 171 of
%   International Series of Numerical Mathematics, Birkhäuser, Cham, 2021,
%   pp. 369–392.
%   DOI: 10.1007/978-3-030-72983-7_18
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
