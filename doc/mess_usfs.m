%
%  ######################################################################
%  #                                                                    #
%  #         M-M.E.S.S. USFS - the user supplied function system        #
%  #                                                                    #
%  ######################################################################
%
%  All operations with the implicit matrices (see (1) below) E, A, and
%  (A + p E) regarding multiplication and solution of linear systems
%  are implemented via user supplied functions used via function
%  handles in the oper structure. The USFS section in `help mess`
%  describes the sets of function handles shipped with M-M.E.S.S. and
%  what types of systems they implement.
%
%  Here, we describe which functions are mandatory, which
%  additional ones are optional and how they can be used with the
%  M-M.E.S.S. 'operatormanager' that generates the corresponding
%  'oper' structure.
%
%  All operations in M-M.E.S.S are implemented with a dynamical
%  system
%      .
%    E x(t) = A x(t) + B u(t),                                 (1)
%      y(t) = C x(t),
%
%  where E is invertible, and all matrices are real, in mind.
%
%  Other system structures can be addressed, as long as they have an
%  equivalent realization of this form that is accessible from the
%  given data. This would, rather obviously, be the case for higher
%  order differential systems, where always companion-form
%  realizations in first order form can be set up as block
%  matrices. For proper differential algebraic cases, (1) represents
%  the implicit projection onto the solution manifold. The key idea
%  behind the usfs system is that forming these matrices can be
%  avoided and solving linear systems, or performing multiplications
%  with these matrices, can be expressed in terms of the original
%  matrices. This is often possible for far more complex
%  situations. Also matrices do not have to be present in the MATLAB
%  workspace, e.g. for calling your favorite discretization tool via
%  the mex interface.
%
%  We require a few mandatory functions for each usfs set with a
%  fixed naming scheme. To this end consider the set is called
%  'setname' residing in the folder 'setname' inside the folder
%  'usfs'
%
%  n = size_setname(eqn, opts, oper)
%     returns the number of rows of A, or its equivalent realization, in (1).
%     for the proper differential algebraic cases this is the dimension of the
%     solution manifold.
%
%  Operations with A:
%
%  C = mul_A_setname(eqn, opts, opA, B, opB)
%     implements the equivalent of
%
%                C = A * B
%                C = A'* B
%                C = A * B'
%             or C = A'* B'
%
%     depending on the values of opA and opB.
%
%  C = sol_A_setname(eqn, opts, opA, B, opB)
%     implements the equivalent of
%
%                C = A \ B
%                C = A'\ B
%                C = A \ B'
%             or C = A'\ B'
%
%     depending on the values of opA and opB.
%
%  Operation with E:
%     analogous to the above, with the specialty,
%     that in the core codes both operations are only called when
%     eqn.haveE is anything that evaluates true.
%
%  Operations with (A + p * E)-type matrices:
%
%  C = mul_ApE_setname(eqn, opts, opA, p, opE, B, opB)
%     implements the equivalent of
%
%                C = (A + p * E ) * B
%                C = (A'+ p * E ) * B
%                C = (A'+ p * E') * B
%                C = (A + p * E') * B
%                C = (A + p * E ) * B'
%                C = (A'+ p * E ) * B'
%                C = (A'+ p * E') * B'
%            or  C = (A + p * E') * B'
%
%     depending on the values of opA, opE and opB.
%
%  C = sol_ApE_setname(eqn, opts, opA, p, opE, B, opB)
%     implements the equivalent of
%
%                C = (A + p * E ) \ B
%                C = (A'+ p * E ) \ B
%                C = (A'+ p * E') \ B
%                C = (A + p * E') \ B
%                C = (A + p * E ) \ B'
%                C = (A'+ p * E ) \ B'
%                C = (A'+ p * E') \ B'
%            or  C = (A + p * E') \ B'
%
%     depending on the values of opA, opE and opB.
%
%  Initialization routines:
%
%  [result, eqn, opts, oper] = init_setname(eqn, opts, oper, flag1, flag2)
%
%     function for the initialization of the set. It essentially
%     checks consistency and availability of all required data for
%     the above operations. Both flag1 and flag2 can be 'E' and
%     'A' to ask for checking the data for the corresponding implicit
%     matrix.
%
%  [RHS, res0, eqn, opts, oper] = init_res_setname(eqn, opts, oper, RHS)
%
%     This function is likely going to get renamed in version 3.0,
%     since it actually serves 2 purposes at the moment. The first
%     of which is exactly what the name says. RHS is the factor of
%     the constant term in the matrix equation on input. res0 is its
%     norm as requested by opts.norm. RHS on output can be the
%     same as the input, but may also be the projection to the
%     coordinates of the implicit representation, i.e. for the
%     differential algebraic cases (the DAE usfs sets) projection
%     to the solution manifold, where the implicit ODE
%     representation (1) evolves. The later motivates the second
%     use of the function as the main tool for explicit projection
%     of data to those implicit representations.
%
%
%  Additional optional functions:
%
%  With the exception of init_setname all usfs have two optional
%  companion functions with suffixes _pre and _post, as in
%
%     [eqn, opts, oper] = sol_ApE_pre_default(eqn, opts, oper)
%
%  These functions are called in each code that relies on the
%  operation. While _pre is used for initialization of data used in
%  the operation (e.g. sol_A_pre could precompute a matrix
%  factorization or a preconditioner) at the very beginning, _post
%  is intended for cleaning up these additional things at the very
%  end of the function.
%
%  [eqn, opts, oper] = eval_matrix_functions_setname(eqn, opts, oper, t)
%
%  In case of time-dependent matrices in (1), i.e. (1) is linear
%  time-varying (LTV) and the corresponding differential matrix
%  equations are non-autonomous, this function evaluates the
%  matrix-valued functions at time t. All following usfs operation use
%  the data at time t until eval_matrix_functions_setname is called
%  again. In the LTV case this function is required and used to
%  initialize the matrices at the start time.
%
%  [rw, Hp, Hm, Vp, Vm] = get_ritz_vals_setname(eqn, opts, oper, U, W,
%                         p_old)
%  This function allows to modify the data before shift parameters are
%  computed, e.g. for the function mess_lradi.m. For the differential
%  algebraic cases the shift parameters are computed with the original data,
%  instead of the transformed data on the solution manifold.
%
%  C = dss_to_ss_state_space_transformed_default...
%      (eqn, opts, fac, opFac, B, opB)
%  C = ss_to_dss_state_space_transformed_default...
%      (eqn, opts, fac, opFac, B, opB)
%
%  This pair of functions is only used in state_space_transformed
%  usfs, where they encode the transformation of (1) to standard
%  statespace form, i.e. an equivalent system with E=I the identity.
%
%  Additional data for the usfs can be stored in eqn, opts or oper, as
%  needed. This is the main reason why all mess_* functions have these
%  as input as well as output arguments. It is usually a good idea to
%  equip the data with counters, such that e.g. matrix factorizations
%  are only computed once in nested function calls, and only get freed
%  once the outermost call wants to delete them. This is done
%  appropriately, where needed, for all usfs that we ship with
%  M-M.E.S.S., already.
%
%  Optional functions that are not present will be linked to
%  mess_do_nothing by the operatormanager.
%
%  Further details on the single sets of usfs can be found via
%
%     help mess_usfs_dae_1
%     help mess_usfs_dae_1_so
%     help mess_usfs_dae_2
%     help mess_usfs_default_iter
%     help mess_usfs_default
%     help mess_usfs_so_1
%     help mess_usfs_so_2
%     help mess_usfs_so_iter
%     help mess_usfs_state_space_transformed_default
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
