function [out, eqn, opts, oper] = mess_lrnm(eqn, opts, oper)
%% function [out, eqn, opts, oper] = mess_lrnm(eqn, opts, oper)
%
% Solve continuous-time Riccati equations with sparse coefficients with
% Newton's method (NM)
%
%   eqn.type = 'N'
%     A*X*E' + E*X*A' - E*X*C'*C*X*E' + B*B' = 0
%   or
%     A*X*E' + E*X*A' - E*X*C'*Q\C*X*E' + B*R*B' = 0
%
%   eqn.type = 'T'
%     A'*X*E + E'*X*A - E'*X*B*B'*X*E + C'*C = 0
%   or
%     A'*X*E + E'*X*A - E'*X*B*R\B'*X*E + C'*Q*C = 0
%
%
% Matrix A can have the form A = Ã + U*V' if U (eqn.U) and V (eqn.V) are
% provided U and V are dense (n x m3) matrices and should satisfy m3 << n
%
%
% The solution is approximated as X = Z*Z', or if opts.LDL_T is true as
% X = L*D*L'
%
% Input/Output
%   eqn         struct contains data for equations
%
%   opts        struct contains parameters for the algorithm
%
%   oper        struct contains function handles for operation
%               with A and E
%
% Output
%   out         struct containing output information
%
% Input fields in struct eqn:
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   eqn.R       dense symmetric and invertible (m1 x m1) matrix
%               (required for LDL^T formulation)
%
%   eqn.Q       dense symmetric (m2 x m2) matrix
%               (required for LDL^T formulation)
%
%   eqn.U       dense (n x m3) matrix U
%               (required if eqn.V is present)
%
%   eqn.V       dense (n x m3) matrix V
%               (required if eqn.U is present)
%
%   eqn.type    possible values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional)
%
%   eqn.haveE   possible values: false, true
%               if haveE == false: matrix E is assumed to be identity
%               (optional)
%
%   eqn.haveUV  possible values: false, true
%               if haveUV = true: U = [U1, U2] and V = [V1, V2]
%               if K or DeltaK are accumulated during the iteration they
%               use only U2 and V2. U1 and V1 can be used for an external
%               rank-k update of the operator.
%               The size of U1 and V1 can be given via eqn.sizeUV1.
%               (optional, default: false)
%
%   eqn.sizeUV1 possible values: nonnegative integer
%               if a stabilizing feedback is given via U = [U1, U2] and
%               V = [V1, V2] in U2 or V2, eqn.widthU1 indicates how
%               many beginning columns of U and V does not be
%               (optional, default: size(eqn.U, 2))
%
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices. For the second order ode types these are given
%   implicitly by the M, D, K matrices stored in eqn.M_, eqn.E_ and eqn.K_,
%   respectively.
%
% Input fields in struct opts:
%   opts.LDL_T                  possible values: false, true
%                               use LDL^T formulation for the RHS and
%                               solution
%                               (optional, default: false)
%
%   opts.norm                   possible values: 2, 'fro'
%                               use 2-norm (2) or Frobenius norm ('fro') to
%                               compute residual and relative change norms
%                               in case projection is used
%                               (opts.nm.projection.freq > 0) norm will
%                               automatically be set to 2
%                               (optional, default: 'fro')
%
%   opts.nm.K0                  possible values: dense 'T': (m1 x n) matrix
%                               or 'N': (m2 x n) matrix
%                               initial stabilizing feedback K
%                               (optional)
%
%   opts.nm.maxiter             possible  values: integer > 0
%                               maximum NM iteration number
%                               (optional, default: 20)
%
%   opts.nm.res_tol             possible values: scalar >= 0
%                               stopping tolerance for the relative NM
%                               residual norm; if res_tol = 0 the relative
%                               residual norm is not evaluated
%                               (optional, default: 0)
%
%   opts.nm.rel_diff_tol        possible values: scalar >= 0
%                               stopping tolerance for the relative
%                               change of the NM solution Z;
%                               if res_tol = 0 the relative
%                               change is not evaluated
%                               (optional, default: 0)
%
%   opts.nm.trunc_tol           possible values: scalar > 0
%                               tolerance for rank truncation of the
%                               low-rank solutions (aka column compression)
%                               (optional, default: eps*n)
%
%   opts.nm.trunc_info          possible values: 0, 1
%                               verbose mode for column compression
%                               (optional, default: 0)
%
%   opts.nm.info                possible  values: 0, 1
%                               turn on (1) or off (0) the status output in
%                               every NM iteration step
%                               (optional, default: 0)
%
%   opts.nm.accumulateRes       possible  values: false, true
%                               accumulate the relative NM residual norm
%                               during the inner ADI iteration
%                               (optional, default: false)
%
%   opts.nm.linesearch          possible  values: false, true
%                               if turned of (false) NM makes full steps; if
%                               turned on (true) a step size 0<=lambda<=2 is
%                               computed
%                               (optional, default: false)
%
%   opts.nm.inexact             possible  values: false, 'linear',
%                               'superlinear', 'quadratic'
%                               the inner ADI uses an adaptive relative ADI
%                               residual norm; with
%                               ||R||: relative NM residual norm,
%                               tau: opts.nm.tau,
%                               j: NM iteration index;
%                               'linear': tau * ||R||,
%                               'superlinear':  ||R|| / (j^3 + 1),
%                               'quadratic': tau / sqrt(||R||)
%                               (optional, default: false)
%
%   opts.nm.store_solfac        possible values: false, true
%                               if turned on (true) the solution factors
%                               computed by the adi are stored in the
%                               out.adi structure
%                               (optional, default false)
%
%   opts.nm.store_debug         possible values: false, true
%                               if turned on (true) the residual factors and
%                               feedback updates from the adi are stored
%                               in the out.adi structure
%                               (optional, default false)
%
%   opts.nm.tau                 possible  values: scalar >= 0
%                               factor for inexact inner ADI iteration
%                               tolerance
%                               (optional, default: 1)
%
%   opts.nm.projection.freq     possible  values: integer >= 0
%                               frequency of the usage of Galerkin
%                               projection acceleration in NM
%                               (optional, default: 0)
%
% For the following fields see mess_galerkin_projection_acceleration:
%
%   opts.nm.projection.ortho
%   opts.nm.projection.meth
%
%   opts.shifts.period          possible  values: integer > 0
%                               number of NM iterations that should pass
%                               until new shifts are computed
%                               (optional, default: 1)
%
% For the following fields see mess_para:
%
%   opts.shifts.num_desired
%   opts.shifts.method
%
% For the following fields see mess_lradi:
%
%   opts.adi.compute_sol_fac
%   opts.adi.accumulateK
%   opts.adi.accumulateDeltaK
%
% If optional input arguments are missing they may be set to default values
% and a 'MESS:control_data' warning is printed. to turn theses warnings off
% use: warning('OFF', 'MESS:control_data')
%
% For LDL^T formulation use opts.LDL_T = true:
%     RHS of Lyapunov Eq. has form W * T * W'
%     Solution Lyapunov Eq. has form L * D * L'
%     L is stored in Z if computed (opts.adi.compute_sol_fac)
%     T (eqn.T in the ADI for the per Newton step Lyapunov equations) is built
%     from Q and R.
%
% Output fields in struct out:
%   out.Z               low-rank solution factor
%
%   out.adi             struct with the output of the all ADI iterations
%
%   out.niter           number of NM iterations
%
%   out.K               feedback matrix
%                       (T): K = B' ZZ' E
%                       (N): K = C ZZ' E
%                       or
%                       (T): K = R \ B' ZDZ' E
%                       (N): K = Q \ C ZDZ' E
%
%   out.D               solution factor for LDL^T formulation
%                       (opts.LDL_T = true)
%
%   out.res             array of relative NM residual norms
%
%   out.rc              array of relative NM change norms
%
%   out.res0            norm of the normalization residual term
%
%   uses operator functions init and mul_E, mul_E_pre, mul_E_post directly
%   and indirectly:
%     size, init, init_res, sol_ApE and mul_E in mess_lradi;
%     size in mess_para
%     size, sol_A, mul_A, sol_E, mul_E in mess_arn
%     mul_A, mul_E in mess_galerkin_projection_acceleration
%     mul_A, mul_E in riccati
%
%   See also mess_lradi, mess_para, mess_galerkin_projection_acceleration,
%   operatormanager.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for ADI Control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isfield(opts, 'adi')) || not(isstruct(opts.adi))
    mess_err(opts, 'control_data',  ...
             'ADI control structure opts.ADI missing.');
end % Single fields are checked below or inside mess_lradi

if not(isfield(opts.adi, 'compute_sol_fac'))
    opts.adi.compute_sol_fac = true;
end
if not(isfield(opts.adi, 'accumulateK'))
    opts.adi.accumulateK = false;
end
if not(isfield(opts.adi, 'accumulateDeltaK'))
    opts.adi.accumulateDeltaK = false;
end

if not(opts.adi.compute_sol_fac || opts.adi.accumulateK || ...
       opts.adi.accumulateDeltaK)
    mess_warn(opts, 'control_data', ...
              ['Either opts.adi.accumulateK or opts.adi.compute_sol_fac or ', ...
               'opts.adi.accumulateDeltaK must be true. Switching to default', ...
               'opts.adi.accumulateDeltaK = true']);
    opts.adi.accumulateDeltaK = true;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');

if not(result)
    mess_err(opts, 'control_data', ...
             'system data is not completely defined or corrupted');
end

[eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);

[eqn, opts, oper] = oper.init_res_pre(eqn, opts, oper);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for shift parameter structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts, 'shifts')) || not(isstruct(opts.shifts))
    mess_warn(opts, 'control_data', ...
              ['shift parameter control structure missing.', ...
               'Switching to default num_desired = 25, num_Ritz = 50, ' ...
               'num_hRitz = 25.']);
    opts.shifts.num_desired = 25;
    opts.shifts.num_Ritz = 50;
    opts.shifts.num_hRitz = 25;
    opts.shifts.method = 'heur';
    opts.shifts.period = 1;

else
    if not(isfield(opts.shifts, 'num_desired')) || ...
       not(isnumeric(opts.shifts.num_desired))
        mess_warn(opts, 'control_data', ...
                  ['Missing or Corrupted opts.shifts.num_desired field.', ...
                   'Switching to default: 25']);
        opts.shifts.num_desired = 25;
    end

    if not(isfield(opts.shifts, 'method'))
        opts.shifts.method = 'heur';
    end

    if (strcmp(opts.shifts.method, 'heur') || ...
        strcmp(opts.shifts.method, 'wachspress')) && ...
       (not(isfield(opts.shifts, 'num_Ritz')) || ...
        not(isnumeric(opts.shifts.num_Ritz)))

        mess_warn(opts, 'control_data', ...
                  ['Missing or Corrupted opts.shifts.num_Ritz field.', ...
                   'Switching to default: 50']);
        opts.shifts.num_Ritz = 50;
    end

    if (strcmp(opts.shifts.method, 'heur') || ...
        strcmp(opts.shifts.method, 'wachspress')) && ...
       (not(isfield(opts.shifts, 'num_hRitz')) || ...
        not(isnumeric(opts.shifts.num_hRitz)))

        mess_warn(opts, 'control_data', ...
                  ['Missing or Corrupted opts.shifts.num_hRitz field.', ...
                   'Switching to default: 25']);
        opts.shifts.num_hRitz = 25;
    end

    if not(isfield(opts.shifts, 'period'))
        opts.shifts.period = 1;
    end

    if not(isfield(opts.shifts, 'wachspress'))
        opts.shifts.wachspress = 'T';
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for Newton control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts, 'nm')) || not(isstruct(opts.nm))
    mess_err(opts, 'control_data', ...
             'Newton control structure opts.nm missing.');
else
    if not(isfield(opts.nm, 'maxiter')) || not(isnumeric(opts.nm.maxiter))
        mess_warn(opts, 'control_data', ...
                  ['Missing or corrupted ''maxiter'' field. ', ...
                   'Switching to default.']);
        opts.nm.maxiter = 20;
    end

    if not(isfield(opts.nm, 'rel_diff_tol')) || ...
            not(isnumeric(opts.nm.rel_diff_tol))
        mess_warn(opts, 'control_data', ...
                  ['Missing or corrupted ''rel_diff_tol'' field. ', ...
                   'Switching to default.']);
        opts.nm.rel_diff_tol = 0;
    end

    if not(isfield(opts.nm, 'res_tol')) || not(isnumeric(opts.nm.res_tol))
        mess_warn(opts, 'control_data', ...
                  ['Missing or corrupted ''res_tol'' field. ', ...
                   'Switching to default.']);
        opts.nm.res_tol = false;
    end

    if not(isfield(opts.nm, 'accumulateRes')) || ...
            not(islogical(opts.nm.accumulateRes))
        mess_warn(opts, 'control_data', ...
                  ['Missing or corrupted ''accumulateRes'' field. ', ...
                   'Switching to default.']);
        opts.nm.accumulateRes = false;
    end

    if opts.nm.accumulateRes
        % need DeltaK
        opts.adi.accumulateDeltaK = true;
    end

    if not(isfield(opts.nm, 'linesearch')) || ...
            not(islogical(opts.nm.linesearch))
        mess_warn(opts, 'control_data', ...
                  ['Missing or corrupted ''linesearch'' field. ', ...
                   'Switching to default.']);
        opts.nm.linesearch = false;
    end

    if not(isfield(opts.nm, 'store_solfac')) || ...
            not(islogical(opts.nm.store_solfac))
        opts.nm.store_solfac = false;
    end

    if not(isfield(opts.nm, 'store_debug')) || ...
       not(islogical(opts.nm.store_debug))
        opts.nm.store_debug = false;
    end

    if opts.nm.linesearch
        % need DeltaK
        opts.adi.accumulateDeltaK = true;
        % need res_tol
        if not(opts.nm.res_tol)
            opts.nm.res_tol = 1e-16;
        end
        alpha = 1e-4;
    end

    if not(isfield(opts.nm, 'inexact'))
        opts.nm.inexact = false;
    end

    if opts.nm.inexact
        if not(isfield(opts.nm, 'tau'))
            opts.nm.tau = 1;
        end
        if not(isfield(opts.adi, 'res_tol'))
            opts.adi.res_tol = 1e-16;
        end
        opts.nm.accumulateRes = true;
        opts.adi.inexact = true;
        opts.adi.accumulateDeltaK = true;
    else
        opts.adi.inexact = false;
    end

    if not(isfield(opts, 'norm')) || ...
       (not(strcmp(opts.norm, 'fro')) && ...
        (not(isnumeric(opts.norm)) ||  not(opts.norm == 2)))

        mess_warn(opts, 'control_data', ...
                  ['Missing or Corrupted norm field. ', ...
                   'Switching to default.']);
        opts.norm = 'fro';
    end

    if not(isfield(opts.nm, 'info')) || not(isnumeric(opts.nm.res_tol))
        mess_warn(opts, 'control_data', ...
                  ['Missing or Corrupted info field. ', ...
                   'Switching to default.']);
        opts.nm.info = 0;
    end
end

if not(isfield(opts.nm, 'trunc_tol'))
    opts.nm.trunc_tol = eps * oper.size(eqn, opts);
end

if not(isfield(opts.nm, 'trunc_info'))
    opts.nm.trunc_info = 0;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for projection triggers and residual control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project = isfield(opts.nm, 'projection') && ...
          isfield(opts.nm.projection, 'freq') && ...
          isnumeric(opts.nm.projection.freq) && ...
          opts.nm.projection.freq;

if project
    opts.adi.compute_sol_fac = true;
    opts.norm = 2;
end

% we need to use iterative residual computation. Let's add the
% corresponding control structure in case it does not already exist.
if project && not(isfield(opts.nm, 'res'))
    mess_warn(opts, 'control_data', ...
              ['Found empty residual control parameters. ', ...
               'Falling back to defaults.']);
    opts.nm.res = [];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isfield(eqn, 'B')) || not(isnumeric(eqn.B))
    mess_err(opts, 'control_data', 'eqn.B is not defined or corrupted');
end

if not(isfield(eqn, 'C')) || not(isnumeric(eqn.C))
    mess_err(opts, 'control_data', 'eqn.C is not defined or corrupted');
end

% make sure the first right hand side is dense so that the resulting factor
% is densely stored.
if issparse(eqn.C)
    eqn.C = full(eqn.C);
end
if issparse(eqn.B)
    eqn.B = full(eqn.B);
end

if not(isfield(eqn, 'type'))
    eqn.type = 'N';
    mess_warn(opts, 'control_data', ...
              ['Unable to determine type of equation. '...
               'Falling back to type ''N''']);
elseif not(eqn.type == 'N') && not(eqn.type == 'T')
    mess_err(opts, 'equation_type', ...
             'Equation type must be either ''T'' or ''N''');
end

if eqn.type == 'T'
    p = size(eqn.C, 1); % number of outputs
    m = size(eqn.B, 2); % number of inputs
else
    m = size(eqn.C, 1); % number of outputs
    p = size(eqn.B, 2); % number of inputs
end

% check whether LDL^T formulation should be used
if not(isfield(opts, 'LDL_T'))
    opts.LDL_T = false;
end

% check for or set proper right hand side
if opts.LDL_T
    % RHS of inner Lyapunov Eq. has form G * S * G'
    % Solution Lyapunov Eq. has form L * D * L'
    % with D Kronecker product of adiout.D and S
    % D is not computed explicitly
    % L is stored in Z if computed (opts.adi.compute_sol_fac)
    % S (eqn.Q or inv(eqn.R)) need to be given
    if not(isfield(eqn, 'Q')) || not(isnumeric(eqn.Q)) || ...
            not(isfield(eqn, 'R')) || not(isnumeric(eqn.R))
        mess_err(opts, 'control_data', ...
                 'eqn.Q or eqn.R is not defined or corrupted');
    end

    if isfield(opts, 'bdf') && isstruct(opts.bdf)
        if not(isfield(opts.bdf, 'tau')) || not(isnumeric(opts.bdf.tau))
            mess_err(opts, 'control_data', ...
                     'opts.bdf.tau is not defined or corrupted');
        end

        if not(isfield(opts.bdf, 'beta')) || not(isnumeric(opts.bdf.beta))
            mess_err(opts, 'control_data', ...
                     'opts.bdf.beta is not defined or corrupted');
        end

        tau_beta = opts.bdf.tau * opts.bdf.beta;
        bdf = true;
    else
        bdf = false;
        tau_beta = 1;
    end
else
    bdf = false;
    if not(isfield(eqn, 'Q')) || not(isnumeric(eqn.Q))
        eqn.Q = eye(size(eqn.C, 1));
    end
    if not(isfield(eqn, 'R')) || not(isnumeric(eqn.R))
        eqn.R = eye(size(eqn.B, 2));
    end
end
%%
% in order not to overwrite important equation data we send a copy to
% the inner iteration with the proper augmentation
adi_eqn = eqn;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rank-k update system data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subroutines we call rely on haveUV being set, so we add it to both
if not(isfield(eqn, 'haveUV')) || isempty(eqn.haveUV) || ...
        not(eqn.haveUV)
    eqn.haveUV = false;
    eqn.sizeUV1 = 0;
    eqn.U       = [];
    eqn.V       = [];
end

if not(isfield(adi_eqn, 'haveUV')) || isempty(adi_eqn.haveUV) || ...
        not(adi_eqn.haveUV)
    adi_eqn.haveUV = false;
    adi_eqn.sizeUV1 = 0;
    adi_eqn.U       = [];
    adi_eqn.V       = [];
else
    if opts.LDL_T
        mess_err(opts, 'control_data', ...
                 ['LDL_T formulation is not compatible with ' ...
                  'external eqn.haveUV option.']);
    end

    if isnumeric(adi_eqn.U) && isnumeric(adi_eqn.V) && ...
       (size(adi_eqn.U, 1) == size(adi_eqn.V, 1)) && ...
       (size(adi_eqn.U, 2) == size(adi_eqn.V, 2))

        if issparse(adi_eqn.V)
            adi_eqn.V = full(adi_eqn.V);
        end
        if issparse(adi_eqn.U)
            adi_eqn.U = full(adi_eqn.U);
        end
    else
        mess_err(opts, 'control_data', ...
                 ['Inappropriate data of low-rank updated operator', ...
                  '(eqn.U and eqn.V)']);
    end
end

% Check for size of constant term in U and V.
if adi_eqn.haveUV
    if not(isfield(adi_eqn, 'sizeUV1')) || isempty(adi_eqn.sizeUV1)
        adi_eqn.sizeUV1 = size(adi_eqn.U, 2);
    else
        mess_assert(opts, isnumeric(adi_eqn.sizeUV1) && ...
                    (adi_eqn.sizeUV1 <= size(adi_eqn.U, 2)), ...
                    'control_data', ...
                    ['Inappropriate size of low-rank updated operator ', ...
                     '(eqn.U and eqn.V)']);
    end
end

% The actual operations in the Lyapunov solver have two low-rank updates:
%
%   * the one coming in from A = F + U * V' here
%   * the one introduced by the Kleinman step using B (or C) and the
%     feedback matrix
%
% We put both into extended U and V in the inner Lyapunov solver and for easier
% reading we introduce index vectors for those
UV_cols = 1:adi_eqn.sizeUV1;
feedb_cols = adi_eqn.sizeUV1 + 1:adi_eqn.sizeUV1 + m;

% Initialize storage for the computed feedback.
if eqn.type == 'T'
    adi_eqn.U = [adi_eqn.U(:, UV_cols), -eqn.B];
    adi_eqn.V = [adi_eqn.V(:, UV_cols), zeros(size(eqn.B))];
else
    adi_eqn.U = [adi_eqn.U(:, UV_cols), zeros(size(eqn.C'))];
    adi_eqn.V = [adi_eqn.V(:, UV_cols), -eqn.C'];
end

adi_eqn.haveUV = true;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if we have an initial stabilizing feedback.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(opts.nm, 'K0')
    if eqn.type == 'T'
        adi_eqn.V(:, feedb_cols) = opts.nm.K0';
    else
        adi_eqn.U(:, feedb_cols) = opts.nm.K0';
    end
end

if bdf
    if eqn.type == 'T'
        adi_eqn.U = -tau_beta * eqn.B;
    else
        adi_eqn.V = -tau_beta * eqn.C';
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All checks done. Here comes the real work!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if changes are made here these have to be repeated at the restart part
% from line search further below!!!
if opts.nm.res_tol
    res = zeros(opts.nm.maxiter, 1);
else
    res = [];
end

if opts.nm.rel_diff_tol
    rc = zeros(opts.nm.maxiter, 1);
else
    rc = [];
end

if eqn.type == 'T'

    adi_eqn.W = eqn.C';
    if opts.LDL_T
        adi_eqn.T = adi_eqn.Q;
    end

else

    adi_eqn.W = eqn.B;
    if opts.LDL_T
        adi_eqn.T = adi_eqn.R;
    end

end

% init_res for projection of data. Actual residual computations follows below.
adi_eqn.W = oper.init_res(adi_eqn, opts, oper, adi_eqn.W);

if opts.LDL_T

    res0 = riccati_LR(adi_eqn.W, [], opts, adi_eqn.T, []);

    if eqn.type == 'T'

        adi_eqn.T = blkdiag(eqn.Q, tau_beta * eqn.R);

    else

        adi_eqn.T = blkdiag(eqn.R, tau_beta * eqn.Q);

    end

else

    res0 = norm(adi_eqn.W' * adi_eqn.W, opts.norm);
    adi_eqn.T = [];

end

opts.nm.res0 = res0;

if eqn.type == 'T'
    adi_eqn.W = [adi_eqn.W, adi_eqn.V(:, feedb_cols)];
else
    adi_eqn.W = [adi_eqn.W, adi_eqn.U(:, feedb_cols)];
end

if opts.adi.inexact
    switch opts.nm.inexact

        case 'linear'
            opts.adi.outer_tol = opts.nm.tau * res0;

        case 'superlinear'
            opts.adi.outer_tol = 0.5 * res0;

        case 'quadratic'
            if res0 > 1
                opts.adi.outer_tol = opts.nm.tau / sqrt(res0);
            else
                opts.adi.outer_tol = opts.nm.tau * res0 * res0;
            end

        otherwise
            mess_err(opts, 'inexact', ...
                     ['inexact must be false, ''linear'', ''superlinear''', ...
                      ' or ''quadratic''']);
    end

    opts.adi.outer_tol = max(opts.adi.outer_tol, opts.adi.res_tol);
end

if opts.nm.linesearch
    linesearch = false;
    W_old = adi_eqn.W;
    DeltaK_old = [];
    if opts.LDL_T
        S_old = adi_eqn.T;
    end
end

restarted = false;

already_restarted = false;
% if changes are made here these have to be repeated at the restart part
% from line search further below!!!

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;

while k <= opts.nm.maxiter

    projected = false;

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute new ADI shifts every "period" steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if not(mod(k - 1, opts.shifts.period))
        opts.shifts.p = mess_para(adi_eqn, opts, oper);
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % form right hand side factor appending the current feedback
    % approximant to the factor of the constant term
    % first p columns of G are always the same and were initialized above
    if eqn.type == 'T'
        adi_eqn.W(:, p + 1:p + m) = adi_eqn.V(:, feedb_cols);
    else
        adi_eqn.W(:, p + 1:p + m) = adi_eqn.U(:, feedb_cols);
    end

    if opts.LDL_T
        if eqn.type == 'T'
            adi_eqn.T = blkdiag(eqn.Q, tau_beta * eqn.R);
        else
            adi_eqn.T = blkdiag(eqn.R, tau_beta * eqn.Q);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve the Lyapunov equation

    [adiout, adi_eqn, opts, oper] = mess_lradi(adi_eqn, opts, oper);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store only necessary information about ADI iteration.
    tmp = adiout;

    if not(opts.nm.store_solfac)
        if isfield(tmp, 'Z')
            tmp = rmfield(tmp, 'Z');
        end
        if isfield(tmp, 'D')
            tmp = rmfield(tmp, 'D');
        end

    end

    if not(opts.nm.store_debug)
        if isfield(tmp, 'DeltaK')
            tmp = rmfield(tmp, 'DeltaK');
        end
        if isfield(tmp, 'Knew')
            tmp = rmfield(tmp, 'Knew');
        end
        if isfield(tmp, 'res_fact')
            tmp = rmfield(tmp, 'res_fact');
        end
    end
    out.adi(k) = tmp;

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Restart Newton iteration with exact ADI iteration if necessary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if adiout.restart && opts.nm.inexact(1)

        if already_restarted
            mess_err(opts, 'lrnm', ...
                     ['Newton iteration with line search ', ...
                      'failed to converge.']);
        else
            % restart Newton iteration
            mess_warn(opts, 'lrnm', ...
                      ['Newton iteration needs to be restarted ', ...
                       'because of divergence. ', ...
                       'Continuing with exact ADI iteration.']);
            if opts.nm.res_tol
                res = [res(1:k - 1); zeros(opts.nm.maxiter, 1)];
            else
                res = [];
            end

            if opts.nm.rel_diff_tol
                rc = [rc(1:k - 1); zeros(opts.nm.maxiter, 1)];
            else
                rc = [];
            end

            opts.nm.maxiter = opts.nm.maxiter + k;

            % Reset U and V to initial state.
            if eqn.type == 'T'
                adi_eqn.U = [adi_eqn.U(:, UV_cols), -adi_eqn.B];
                adi_eqn.V = [adi_eqn.V(:, UV_cols), zeros(size(adi_eqn.B))];
            else
                adi_eqn.U = [adi_eqn.U(:, UV_cols), zeros(size(adi_eqn.C'))];
                adi_eqn.V = [eqn.V(:, UV_cols), -adi_eqn.C'];
            end

            adi_eqn.haveUV = true;

            if isfield(opts.nm, 'K0')
                if eqn.type == 'T'
                    adi_eqn.V(:, feedb_cols) = opts.nm.K0';
                else
                    adi_eqn.U(:, feedb_cols) = opts.nm.K0';
                end
            end

            if bdf
                if eqn.type == 'T'
                    adi_eqn.U = -tau_beta * eqn.B;
                else
                    adi_eqn.V = -tau_beta * eqn.C';
                end
            end

            if eqn.type == 'T'

                adi_eqn.W = adi_eqn.C';
                if opts.LDL_T
                    adi_eqn.T = adi_eqn.Q;
                end
            else

                adi_eqn.W = adi_eqn.B;
                if opts.LDL_T
                    adi_eqn.T = adi_eqn.R;
                end

            end

            % use init res for projection, no T required also in LDL_T mode
            adi_eqn.W = oper.init_res(adi_eqn, opts, oper, adi_eqn.W);

            if opts.LDL_T

                if eqn.type == 'T'
                    res0 = riccati_LR(adi_eqn.W, [], opts, adi_eqn.Q, []);
                    adi_eqn.T = blkdiag(eqn.Q, tau_beta * eqn.R);
                else
                    res0 = riccati_LR(adi_eqn.W, [], opts, adi_eqn.R, []);
                    adi_eqn.T = blkdiag(eqn.R, tau_beta * eqn.Q);
                end

            else

                res0 = norm(adi_eqn.W' * adi_eqn.W, opts.norm);

            end

            opts.nm.res0     = res0;
            opts.nm.inexact  = false;
            opts.adi.inexact = false;

            if eqn.type == 'T'
                adi_eqn.W = [adi_eqn.W, adi_eqn.V(:, feedb_cols)];
            else
                adi_eqn.W = [adi_eqn.W, adi_eqn.U(:, feedb_cols)];
            end

            if opts.nm.linesearch
                linesearch = false;
                W_old = adi_eqn.W;
                DeltaK_old = [];
                if opts.LDL_T
                    S_old = adi_eqn.T;
                end
            end
            restarted = true;
            continue
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform projection update if desired
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if project && not(mod(k, opts.nm.projection.freq)) && ...
       not(opts.nm.accumulateRes && ...
           (adiout.Riccati_res < opts.nm.res_tol))

        projected = true;

        opts.nm.linesearch = false; % no line search after projection

        opts.nm.projection.Z = adiout.Z;
        if opts.LDL_T
            opts.nm.projection.D = adiout.D;
        end

        [eqn, opts, oper] = ...
            mess_solve_projected_eqn(eqn, opts, oper, 'GPA', 'CARE');

        adiout.Z = opts.nm.projection.Z;

        if opts.LDL_T
            adiout.D = opts.nm.projection.D;
        end

    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update feedback approximation if it has not been accumulated during
    % the inner iteration or after projection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if not(opts.adi.accumulateK) || projected

        if not(opts.adi.accumulateDeltaK) || projected

            [adiout, eqn, opts, oper] = ...
                mess_accumulateK(eqn, opts, oper, adiout, [], adiout.Z);
        else
            if eqn.type == 'T'
                adiout.Knew = adi_eqn.V(:, feedb_cols) + adiout.DeltaK;
            else
                adiout.Knew = adi_eqn.U(:, feedb_cols) + adiout.DeltaK;
            end
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.nm.res_tol || opts.nm.rel_diff_tol
        if eqn.type == 'T'
            V1 = adiout.Knew - adi_eqn.V(:, feedb_cols);
        else
            V1 = adiout.Knew - adi_eqn.U(:, feedb_cols);
        end
    end

    if opts.nm.res_tol
        if projected
            if opts.LDL_T
                res(k) = mess_res2_norms(adiout.Z, 'riccati', ...
                                         eqn, opts, oper, ...
                                         opts.nm, adiout.D) / res0;
            else
                res(k) = mess_res2_norms(adiout.Z, 'riccati', ...
                                         eqn, opts, oper, ...
                                         opts.nm, []) / res0;
            end
        else
            if opts.nm.accumulateRes
                res(k) = adiout.Riccati_res;
            else
                res(k) = riccati_LR(adiout.res_fact, V1, opts, ...
                                    adi_eqn.T, []) / res0;
            end
        end
    end

    if opts.nm.linesearch
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check whether line search is necessary
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ((k == 1) && (res(k) > 1)) || ...
                ((k > 1) && (res(k) > res(k - 1))) || ...
                adiout.linesearch

            linesearch = true;

            % Compute Lambda
            if opts.LDL_T
                if eqn.type == 'T'
                    M = eqn.R;
                else
                    M = eqn.Q;
                end
                lambda = exact_line_search(opts, ...
                                           W_old, DeltaK_old, ...
                                           adiout.res_fact, ...
                                           adiout.DeltaK, ...
                                           adi_eqn.T, ...
                                           S_old, ...
                                           M);
            else
                lambda = exact_line_search(opts, W_old, DeltaK_old, ...
                                           adiout.res_fact, ...
                                           adiout.DeltaK, [], [], []);
            end

            if opts.nm.info
                mess_fprintf(opts, ...
                             ['\n\t\t Using line search (res: %4d)\n', ...
                              '\t\t lambda: %e\n'], res(k), lambda);
            end

            % Update K
            if eqn.type == 'T'
                adiout.Knew = adi_eqn.V(:, feedb_cols) + ...
                    lambda * adiout.DeltaK;
            else
                adiout.Knew = adi_eqn.U(:, feedb_cols) + ...
                    lambda * adiout.DeltaK;
            end

            % Update DeltaK and W
            if lambda <= 1
                adiout.DeltaK = [sqrt(1 - lambda) * DeltaK_old, ...
                                 lambda * adiout.DeltaK];
                adiout.res_fact = [sqrt(1 - lambda) * W_old, ...
                                   sqrt(lambda) * adiout.res_fact];
                if opts.LDL_T
                    S_linesearch = blkdiag(S_old, adi_eqn.T);
                    S_K = [];
                end
            else
                W_new = adiout.res_fact;
                DeltaK_new = adiout.DeltaK;
                sqrt_lambda = sqrt(lambda);
                adiout.DeltaK = [W_old, sqrt_lambda * W_new, ...
                                 sqrt_lambda * DeltaK_old];
                adiout.res_fact = [DeltaK_old, sqrt_lambda * W_old, ...
                                   lambda * DeltaK_new];

                if opts.LDL_T
                    if eqn.type == 'N'
                        D_K = eqn.Q;
                    else
                        D_K = eqn.R;
                    end

                    if bdf
                        D_K = tau_beta * D_K;
                    end

                    S_linesearch = blkdiag(D_K, S_old, D_K);
                    S_K = blkdiag(S_old, adi_eqn.T, D_K);

                end
            end

            if not(opts.LDL_T)
                S_linesearch = adi_eqn.T;
                S_K = [];
            end

            % Compute residual norm after line search
            res(k) = riccati_LR(adiout.res_fact, adiout.DeltaK, opts, ...
                                S_linesearch, S_K) / res0;
            if k == 1
                bound = (1 - lambda * alpha) * res0;
            else
                bound = (1 - lambda * alpha) * res(k - 1);
            end

            if not(restarted) && (res(k) >= bound)

                % No sufficient decrease
                mess_warn(opts, 'lrnm', ...
                          ['Newton iteration with line search has', ...
                           ' insufficient decrease in iteration k = %d\n', ...
                           '(%g >= %g, LDL_T = %d, eqn.type = %s ', ...
                           'inexact ADI = %s shifts = %s).'], ...
                          k, res(k), bound, opts.LDL_T, eqn.type, ...
                          mess_string(opts.adi.inexact), opts.shifts.method);

                if opts.adi.inexact
                    % switch to exact ADI iteration
                    mess_warn(opts, 'lrnm', ...
                              'Switching to exact ADI iteration.');
                    opts.adi.inexact = false;
                else
                    mess_err(opts, 'lrnm', ...
                             ['Newton iteration with line search ' ...
                              'failed to converge.']);
                end
            end
        end

        % Keep ADI LR residual factor and DeltaK for next Newton iteration
        W_old = adiout.res_fact;
        DeltaK_old = adiout.DeltaK;

        if opts.LDL_T
            if linesearch
                S_old = S_linesearch;
            else
                S_old = adi_eqn.T;
            end
        end

        linesearch = false;
    end

    if opts.nm.rel_diff_tol
        if opts.adi.accumulateDeltaK
            if projected
                if eqn.type == 'T'
                    adiout.DeltaK = adiout.Knew - adi_eqn.V(:, feedb_cols);
                else
                    adiout.DeltaK = adiout.Knew - adi_eqn.U(:, feedb_cols);
                end
            end
            rc(k) = norm(adiout.DeltaK, opts.norm) / ...
                    norm(adiout.Knew, opts.norm);
        else
            rc(k) = norm(V1, opts.norm);
        end
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.nm.info
        if opts.nm.rel_diff_tol && opts.nm.res_tol
            mess_fprintf(opts, ...
                         ['\n NM step: %4d normalized residual: \t%e\n' ...
                          '               relative change in K: \t%e\n' ...
                          '               number of ADI steps: \t%d \n\n'], ...
                         k, res(k), rc(k), adiout.niter);
        elseif opts.nm.res_tol
            mess_fprintf(opts, ...
                         ['\n NM step: %4d normalized residual: \t%e \n\n' ...
                          '                number of ADI steps: \t%d \n\n'], ...
                         k, res(k), adiout.niter);
        elseif opts.nm.rel_diff_tol
            mess_fprintf(opts, ...
                         ['\n NM step: %4d relative change in K: \t%e\n\n' ...
                          '               number of ADI steps: \t%d \n\n'], ...
                         k, rc(k));
        end
    end

    if eqn.type == 'T'
        adi_eqn.V(:, feedb_cols) = adiout.Knew;
    else
        adi_eqn.U(:, feedb_cols) = adiout.Knew;
    end

    k = k + 1;

    if restarted
        already_restarted = true;
        restarted = false;
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (opts.nm.res_tol && (res(k - 1) < opts.nm.res_tol)) || ...
       (opts.nm.rel_diff_tol && (rc(k - 1) < opts.nm.rel_diff_tol))
        break
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set tolerance for next ADI iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if not(opts.adi.inexact) && (k > 2)
        % Check whether inexact ADI iteration can be turned on again
        switch opts.nm.inexact
            case 'linear'
                opts.adi.inexact = res(k - 1) < opts.nm.tau * res(k - 2);
            case 'superlinear'
                opts.adi.inexact = res(k - 1) < 1 / (k^3 + 1) * res(k - 2);
            case 'quadratic'
                if res(k - 2) > 1
                    opts.adi.inexact = ...
                        res(k - 1) < opts.nm.tau / sqrt(res(k - 2));
                else
                    opts.adi.inexact = ...
                        res(k - 1) < opts.nm.tau * res(k - 2) * res(k - 2);
                end
        end

        if opts.adi.inexact
            mess_warn(opts, 'lrnm', ...
                      'Turning inexact ADI iteration back on.');
        end
    end

    if opts.adi.inexact
        switch opts.nm.inexact
            case 'linear'
                opts.adi.outer_tol = opts.nm.tau * res(k - 1);
            case 'superlinear'
                opts.adi.outer_tol = 1 / (k^3 + 1) * res(k - 1);
            case 'quadratic'
                if res(k - 1) > 1
                    opts.adi.outer_tol = opts.nm.tau / sqrt(res(k - 1));
                else
                    opts.adi.outer_tol = opts.nm.tau * res(k - 1)^2;
                end
        end
        opts.adi.outer_tol = max(opts.adi.outer_tol, opts.adi.res_tol);
        opts.adi.outer_res = res(k - 1);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.niter = k - 1;
if eqn.type == 'T'
    out.K = adi_eqn.V(:, feedb_cols)';
else
    out.K = adi_eqn.U(:, feedb_cols)';
end

if opts.LDL_T && isfield(adiout, 'Z') && not(isempty(adiout.Z))
    out.D = adiout.D;
    out.Z = adiout.Z;
else
    if isfield(adiout, 'Z') && not(isempty(adiout.Z))
        out.Z = mess_column_compression(adiout.Z, 'N', [], ...
                                        opts.nm.trunc_tol, ...
                                        opts.nm.trunc_info);
    end
end

if opts.nm.res_tol
    out.res = res(1:out.niter);
end

if opts.nm.rel_diff_tol
    out.rc = rc(1:out.niter);
end

% Delete the res0 part from opts for later use of the struct.
if isfield(opts.nm, 'res0')
    out.res0 = opts.nm.res0;
    opts.nm  = rmfield(opts.nm, 'res0');
elseif exist('res0', 'var')
    out.res0 = res0;
end

if out.niter == opts.nm.maxiter
    mess_warn(opts, 'convergence', ...
              ['LR-NM reached maximum iteration number. ', ...
               'Results may be inaccurate!']);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);

[eqn, opts, oper] = oper.init_res_post(eqn, opts, oper);
