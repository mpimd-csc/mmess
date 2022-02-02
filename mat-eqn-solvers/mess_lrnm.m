function [out, eqn, opts, oper] = mess_lrnm(eqn, opts, oper)
%% function [out, eqn, opts, oper] = mess_lrnm(eqn, opts, oper)
%
% Solve continuous-time Riccati equations with sparse coefficients with
% Newton's method (NM)
%   eqn.type = 'N' -> A*Z*Z'*E' + E*Z*Z'*A' - E*Z*Z'*C'*C*Z*Z'*E' + B*B' = 0 (N)
%   eqn.type = 'T' -> A'*Z*Z'*E + E'*Z*Z'*A - E'*Z*Z'*B*B'*Z*Z'*E + C'*C = 0 (T)
%
%
% Matrix A can have the form A = Ã + U*V' if U (eqn.U) and V (eqn.V) are
% provided U and V are dense (n x m3) matrices and should satisfy m3 << n
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
%   eqn.G       dense (n x m1) matrix G
%               if present it is used instead of B as RHS
%               (required for LDL^T formulation otherwise optional)
%
%   eqn.S       dense (m1 x m1) matrix (N) or (m2 x m2) matrix (T)
%               expected to be symmetric
%               (required for LDL^T formulation)
%
%   eqn.U       dense (n x m3) matrix U
%               (required if eqn.V is present)
%
%   eqn.V       dense (n x m3) matrix V
%               (required if eqn.U is present)
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional)
%
%   eqn.haveE   possible  values: 0, 1, false, true
%               if haveE = 0: matrix E in eqn.E_ is assumed to be identity
%               (optional)
%
%   eqn.haveUV  possible  values: 0, 1, false, true
%               if haveUV = 1: U = [U1, U2] and V = [V1, V2]
%               if K or DeltaK are accumulated during the iteration they
%               use only U2 and V2. U1 and V1 can be used for an external
%               rank-k update of the operator.
%               The size of U1 and V1 can be given via eqn.sizeUV1.
%               (optional, default: 0)
%
%   eqn.sizeUV1 possible values: nonnegative integer
%               if a stabilizing feedback is given via U = [U1, U2] and
%               V = [V1, V2] in U2 or V2, eqn.widthU1 indicates how
%               many beginning columns of U and V does not be
%               (optional, default: size(eqn.U, 2))
%
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices. For the second order types these are given
%   implicitly by the M, D, K matrices stored in eqn.M_, eqn.E_ and eqn.K_,
%   respectively.
%
% Input fields in struct opts:
%   opts.LDL_T                  possible  values: 0, 1, false, true
%                               use LDL^T formulation for the RHS and
%                               solution
%                               (optional, default: 0)
%
%   opts.norm                   possible  values: 2, 'fro'
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
%   opts.nm.res_tol             possible  values: scalar >= 0
%                               stopping tolerance for the relative NM
%                               residual norm; if res_tol = 0 the relative
%                               residual norm is not evaluated
%                               (optional, default: 0)
%
%   opts.nm.rel_diff_tol        possible  values: scalar >= 0
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
%   opts.nm.trunc_info          possible values: 0, 1, false, true
%                               verbose mode for column compression
%                               (optional, default: 0)
%
%   opts.nm.info                possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the status output in
%                               every NM iteration step
%                               (optional, default: 0)
%
%   opts.nm.accumulateRes       possible  values: 0, 1, false, true
%                               accumulate the relative NM residual norm
%                               during the inner ADI iteration
%                               (optional, default: 0)
%
%   opts.nm.linesearch          possible  values: 0, 1, false, true
%                               if tuned of (0) NM makes full steps; if
%                               turned on (1) a step size 0<=lambda<=2 is
%                               computed
%                               (optional, default: 0)
%
%   opts.nm.inexact             possible  values: 0, false, 'linear',
%                               'superlinear', 'quadratic'
%                               the inner ADI uses an adaptive relative ADI
%                               residual norm; with
%                               ||R||: relative NM residual norm,
%                               tau: opts.nm.tau,
%                               j: NM iteration index;
%                               'linear': tau * ||R||,
%                               'superlinear':  ||R|| / (j^3 + 1),
%                               'quadratic': tau / sqrt(||R||)
%                               (optional, default: 0)
%
%   opts.nm.store_solfac        possible values: 0, 1, false, true
%                               if turned on (1) the solution factors
%                               computed by the adi are stored in the
%                               out.adi structure
%                               (optional, default 0)
%
%   opts.nm.store_debug         possible values: 0, 1, false, true
%                               if turned on (1) the residual factors and
%                               feedback updates from the adi are stored
%                               in the out.adi structure
%                               (optional, default 0)
%
%   opts.nm.tau                 possible  values: scalar >= 0
%                               factor for inexact inner ADI iteration
%                               tolerance
%                               (optional, default: 1)
%
%   opts.nm.projection.freq     possible  values: integer >= 0
%                               frequency of the usage of galerkin
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
% For LDL^T formulation use opts.LDL_T = 1:
%     RHS of Lyapunov Eq. has form G * S * G'
%     Solution Lyapunov Eq. has form L * D * L'
%     with D Kronecker product of adiout.D and S
%     L is stored in Z if computed (opts.adi.compute_sol_fac)
%     S (eqn.S) needs to be given
%
% Output fields in struct out:
%   out.Z               low rank solution factor
%
%   out.adi             struct with the output of the all ADI iterations
%
%   out.niter           number of NM iterations
%
%   out.K               feedback matrix 
%                       (T): K = B' ZZ' E 
%                       (N): K = C ZZ' E
%
%   out.D               solution factor for LDL^T formulation
%                       (opts.LDL_T = 1)
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
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for ADI Control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isfield(opts,'adi')) || not(isstruct(opts.adi))
    error('MESS:control_data','ADI control structure opts.ADI missing.');
end % Single fields are checked below or inside mess_lradi

if not(isfield(opts.adi,'compute_sol_fac')), opts.adi.compute_sol_fac = 1; end
if not(isfield(opts.adi,'accumulateK')), opts.adi.accumulateK = 0; end
if not(isfield(opts.adi,'accumulateDeltaK')), opts.adi.accumulateDeltaK = 0; end

if not(opts.adi.compute_sol_fac || opts.adi.accumulateK || opts.adi.accumulateDeltaK)
    warning('MESS:control_data', ...
        ['Either opts.adi.accumulateK or opts.adi.compute_sol_fac or ', ...
        'opts.adi.accumulateDeltaK must be 1. Switching to default', ...
        'opts.adi.accumulateDeltaK = 1']);
    opts.adi.accumulateDeltaK = 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for shift parameter structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts,'shifts')) || not(isstruct(opts.shifts))
    warning('MESS:control_data',...
            ['shift parameter control structure missing.', ...
             'Switching to default num_desired = 25, num_Ritz = 50, ' ...
             'num_hRitz = 25.']);
    opts.shifts.num_desired = 25;
    opts.shifts.num_Ritz = 50;
    opts.shifts.num_hRitz = 25;
    opts.shifts.method='heur';
    opts.shifts.period=1;

else
    if not(isfield(opts.shifts,'num_desired')) ||...
       not(isnumeric(opts.shifts.num_desired))
        warning('MESS:control_data', ...
                ['Missing or Corrupted opts.shifts.num_desired field.', ...
                 'Switching to default: 25']);
        opts.shifts.num_desired = 25;
    end

    if not(isfield(opts.shifts,'method'))
        opts.shifts.method = 'heur';
    end

    if (strcmp(opts.shifts.method,'heur') || ...
        strcmp(opts.shifts.method,'wachspress')) && ...
       (not(isfield(opts.shifts,'num_Ritz')) || ...
        not(isnumeric(opts.shifts.num_Ritz)))

        warning('MESS:control_data', ...
                ['Missing or Corrupted opts.shifts.num_Ritz field.', ...
                'Switching to default: 50']);
        opts.shifts.num_Ritz = 50;
    end

    if (strcmp(opts.shifts.method,'heur') || ...
        strcmp(opts.shifts.method,'wachspress')) && ...
       (not(isfield(opts.shifts,'num_hRitz')) || ...
        not(isnumeric(opts.shifts.num_hRitz)))

        warning('MESS:control_data',...
                ['Missing or Corrupted opts.shifts.num_hRitz field.', ...
                 'Switching to default: 25']);
        opts.shifts.num_hRitz = 25;
    end

    if not(isfield(opts.shifts,'period'))
        opts.shifts.period = 1;
    end

    if not(isfield(opts.shifts,'wachspress'))
        opts.shifts.wachspress = 'T';
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for projection triggers and residual control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project = isfield(opts.nm,'projection') && ...
          isfield(opts.nm.projection,'freq') && ...
          isnumeric(opts.nm.projection.freq) && ...
          opts.nm.projection.freq;

if project
    opts.adi.compute_sol_fac = 1;
    opts.norm = 2;
end

% we need to use iterative residual computation. Let's add the
% corresponding control structure in case it does not already exist.
if project && not(isfield(opts.nm,'res'))
    warning('MESS:control_data', ...
            'Found empty residual control parameters. Falling back to defaults.');
    opts.nm.res = [];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for Newton control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts,'nm')) || not(isstruct(opts.nm))
    error('MESS:control_data','Newton control structure opts.nm missing.');
else
    if not(isfield(opts.nm,'maxiter')) || not(isnumeric(opts.nm.maxiter))
        warning('MESS:control_data', ...
                'Missing or corrupted ''maxiter'' field. Switching to default.');
        opts.nm.maxiter = 20;
    end

    if not(isfield(opts.nm,'rel_diff_tol')) || not(isnumeric(opts.nm.rel_diff_tol))
        warning('MESS:control_data', ...
                'Missing or corrupted ''rel_diff_tol'' field. Switching to default.');
        opts.nm.rel_diff_tol = 0;
    end

    if not(isfield(opts.nm,'res_tol')) || not(isnumeric(opts.nm.res_tol))
        warning('MESS:control_data', ...
                'Missing or corrupted ''res_tol'' field. Switching to default.');
        opts.nm.res_tol = 0;
    end

    if not(isfield(opts.nm,'accumulateRes')) || not(isnumeric(opts.nm.accumulateRes))
        warning('MESS:control_data', ...
                'Missing or corrupted ''accumulateRes'' field. Switching to default.');
        opts.nm.accumulateRes = 0;
    end

    if opts.nm.accumulateRes
        % need DeltaK
        opts.adi.accumulateDeltaK = 1;
    end

    if not(isfield(opts.nm, 'linesearch')) || not(isnumeric(opts.nm.linesearch))
        warning('MESS:control_data', ...
                'Missing or corrupted ''linesearch'' field. Switching to default.');
        opts.nm.linesearch = 0;
    end

    if not(isfield(opts.nm, 'store_solfac')) || not(isnumeric(opts.nm.store_solfac))
        opts.nm.store_solfac = 0;
    end

    if not(isfield(opts.nm, 'store_debug')) || ...
       not(isnumeric(opts.nm.store_debug))
        opts.nm.store_debug = 0;
    end

    if opts.nm.linesearch
        % need DeltaK
        opts.adi.accumulateDeltaK = 1;
        % need res_tol
        if not(opts.nm.res_tol)
            opts.nm.res_tol = 1e-16;
        end
        alpha = 1e-4;
    end

    if not(isfield(opts.nm,'inexact')),  opts.nm.inexact = false; end

    if opts.nm.inexact
        if not(isfield(opts.nm,'tau')), opts.nm.tau = 1; end
        if not(isfield(opts.adi,'res_tol')),  opts.adi.res_tol = 1e-16; end
        opts.nm.accumulateRes = 1;
        opts.adi.inexact = true;
        opts.adi.accumulateDeltaK = 1;
    else
        opts.adi.inexact = false;
    end

    if not(isfield(opts, 'norm')) || ...
       (not(strcmp(opts.norm, 'fro')) && ...
        (not(isnumeric(opts.norm)) ||  opts.norm ~= 2))

        warning('MESS:control_data', ...
                'Missing or Corrupted norm field. Switching to default.');
        opts.norm = 'fro';
    end

    if not(isfield(opts.nm,'info')) || not(isnumeric(opts.nm.res_tol))
        warning('MESS:control_data',...
                'Missing or Corrupted info field. Switching to default.');
        opts.nm.info = 0;
    end
end

if not(isfield(opts.nm,'trunc_tol'))
    opts.nm.trunc_tol = eps * oper.size(eqn, opts);
end

if not(isfield(opts.nm, 'trunc_info')), opts.nm.trunc_info = 0; end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');

if not(result)
    error('MESS:control_data', ...
          'system data is not completely defined or corrupted');
end

[eqn,opts,oper] = oper.mul_E_pre(eqn,opts,oper);

[eqn, opts, oper] = oper.init_res_pre(eqn, opts, oper);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isfield(eqn, 'B')) || not(isnumeric(eqn.B))
    error('MESS:control_data', 'eqn.B is not defined or corrupted');
end

if not(isfield(eqn, 'C')) || not(isnumeric(eqn.C))
    error('MESS:control_data', 'eqn.C is not defined or corrupted');
end

% make sure the first right hand side is dense so that the resulting factor
% is densely stored.
if issparse(eqn.C), eqn.C = full(eqn.C); end
if issparse(eqn.B), eqn.B = full(eqn.B); end

if not(isfield(eqn, 'type'))
    eqn.type = 'N';
    warning('MESS:control_data',['Unable to determine type of equation.'...
            'Falling back to type ''N''']);
elseif (eqn.type ~= 'N') && (eqn.type ~= 'T')
    error('MESS:equation_type', 'Equation type must be either ''T'' or ''N''');
end

if eqn.type == 'T'
    p = size(eqn.C,1); % number of outputs
    m = size(eqn.B,2); % number of inputs
else
    m = size(eqn.C,1); % number of outputs
    p = size(eqn.B,2); % number of inputs
end

% check whether LDL^T formulation should be used
if not(isfield(opts, 'LDL_T')), opts.LDL_T = 0; end

% check for or set proper right hand side
if opts.LDL_T
    % RHS of Lyapunov Eq. has form G * S * G'
    % Solution Lyapunov Eq. has form L * D * L'
    % with D Kronecker product of adiout.D and S
    % D is not computed explicitly
    % L is stored in Z if computed (opts.adi.compute_sol_fac)
    % S (eqn.S) need to be given
    if not(isfield(eqn, 'S')) || not(isnumeric(eqn.S))
        error('MESS:control_data', 'eqn.S is not defined or corrupted');
    end

    if isdiag(eqn.S)
        out.S = eqn.S;
        eqn.S_diag = diag(eqn.S);
        eqn.diagonalized_RHS = 0;
    else
        % diagonalze S and use U to transform the initial RHS later
        [eqn.U_diag, eqn.S_diag] = eig(eqn.S);
        eqn.S_diag = diag(eqn.S_diag);
        eqn.diagonalized_RHS = 1;
    end

    if isfield(opts, 'bdf') && isstruct(opts.bdf)
        if not(isfield(opts.bdf, 'tau')) || not(isnumeric(opts.bdf.tau))
            error('MESS:control_data', 'opts.bdf.tau is not defined or corrupted');
        end

        if not(isfield(opts.bdf, 'beta')) || not(isnumeric(opts.bdf.beta))
            error('MESS:control_data', 'opts.bdf.beta is not defined or corrupted');
        end

        tau_beta = opts.bdf.tau * opts.bdf.beta;
        bdf = 1;
    else
        bdf = 0;
        tau_beta = 1;
    end
else
    bdf = 0;
    eqn.diagonalized_RHS = 0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rank-k update system data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'haveUV')) || isempty(eqn.haveUV) || not(eqn.haveUV)
    eqn.haveUV  = 0;
    eqn.sizeUV1 = 0;
    eqn.U       = [];
    eqn.V       = [];
else
    if opts.LDL_T
        error('MESS:control_data', ...
              ['LDL_T formulation is not compatible with ' ...
               'external eqn.haveUV option.']);
    end

    if isnumeric(eqn.U) && isnumeric(eqn.V) && ...
       (size(eqn.U, 1) == size(eqn.V, 1)) && (size(eqn.U, 2) == size(eqn.V, 2))

        if issparse(eqn.V), eqn.V = full(eqn.V); end
        if issparse(eqn.U), eqn.U = full(eqn.U); end
    else
        error('MESS:control_data', ...
              'Inappropriate data of low rank updated operator (eqn.U and eqn.V)');
    end
end

% Check for size of constant term in U and V.
if eqn.haveUV
    if not(isfield(eqn, 'sizeUV1')) || isempty(eqn.sizeUV1)
        eqn.sizeUV1 = size(eqn.U, 2);
    else
        assert(isnumeric(eqn.sizeUV1) && (eqn.sizeUV1 <= size(eqn.U, 2)), ...
               'MESS:control_data', ...
               'Inappropriate size of low rank updated operator (eqn.U and eqn.V)');
    end
end

% Initialize storage for the computed feedback.
if eqn.type == 'T'
    eqn.U = [eqn.U(:, 1:eqn.sizeUV1), -eqn.B];
    eqn.V = [eqn.V(:, 1:eqn.sizeUV1), zeros(size(eqn.B))];
else
    eqn.U = [eqn.U(:, 1:eqn.sizeUV1), zeros(size(eqn.C,2), size(eqn.C,1))];
    eqn.V = [eqn.V(:, 1:eqn.sizeUV1), -eqn.C'];
end

eqn.haveUV = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if we have an initial stabilizing feedback.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(opts.nm,'K0')
    if eqn.type == 'T'
        eqn.V(:, eqn.sizeUV1+1:eqn.sizeUV1+m) = opts.nm.K0';
    else
        eqn.U(:, eqn.sizeUV1+1:eqn.sizeUV1+m) = opts.nm.K0';
    end
end

if bdf
    if eqn.type == 'T'
        eqn.U = (-opts.bdf.tau * opts.bdf.beta) * eqn.B;
    else
        eqn.V = (-opts.bdf.tau * opts.bdf.beta) * eqn.C';
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
    res = zeros(opts.nm.maxiter,1);
else
    res = [];
end

if opts.nm.rel_diff_tol
    rc = zeros(opts.nm.maxiter,1);
else
    rc = [];
end


% in the LDL_T case we may have diagonalized the kernel matrix of the RHS.
% If so, we need to initialize the residual with the updated G matrix
if eqn.type == 'T'
    if eqn.diagonalized_RHS
        eqn.G = eqn.C' * eqn.U_diag;
    else
        eqn.G = eqn.C';
    end
else
    if eqn.diagonalized_RHS
        eqn.G = eqn.B * eqn.U_diag;
    else
        eqn.G = eqn.B;
    end
end

eqn.G = oper.init_res(eqn, opts, oper, eqn.G);

if opts.LDL_T
    res0 = riccati_LR(eqn.G, [], opts, diag(eqn.S_diag), []);
    S = eqn.S_diag;
    if eqn.type == 'T'
        eqn.S_diag = [eqn.S_diag; tau_beta * ones(size(eqn.V, 2), 1)];
    else
        eqn.S_diag = [eqn.S_diag; tau_beta * ones(size(eqn.U, 2), 1)];
    end
else
    res0 = norm(eqn.G'*eqn.G, opts.norm);
    eqn.S_diag = [];
end

opts.nm.res0 = res0;

if eqn.type == 'T'
    eqn.G = [eqn.G, eqn.V(:, eqn.sizeUV1+1:end)];
else
    eqn.G = [eqn.G, eqn.U(:, eqn.sizeUV1+1:end)];
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
            error('MESS:inexact', ...
                  ['inexact must be 0, ''linear'', ''superlinear''', ...
                   ' or ''quadratic''']);
    end

    opts.adi.outer_tol = max(opts.adi.outer_tol, opts.adi.res_tol);
end

if opts.nm.linesearch
    linesearch = 0;
    W_old = eqn.G;
%    if not(isfield(eqn,'K0')) % FIXME: dead code
        DeltaK_old = [];
%     else
%         %         DeltaK_old = -eqn.V( : , s + 1 : end);
%         %         careful with haveUV option!
%         DeltaK_old = [];
%    end
    if opts.LDL_T
        S_old = eqn.S_diag;
    end
end

restarted = 0;

already_restarted = 0;
% if changes are made here these have to be repeated at the restart part
% from line search further below!!!

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;

while k <= opts.nm.maxiter

    projected = 0;

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute new ADI shifts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if not(mod(k - 1,opts.shifts.period))
        opts.shifts.p = mess_para(eqn,opts,oper);
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % form right hand side factor
    if eqn.type == 'T'
        eqn.G(:, p+1:p+m) = eqn.V(:, end-m+1:end);
    else
        eqn.G(:, p+1:p+m) = eqn.U(:, end-m+1:end);
    end

    if opts.LDL_T
        if eqn.type == 'T'
            eqn.S_diag = [S; tau_beta * ones(size(eqn.V, 2), 1)];
        else
            eqn.S_diag = [S; tau_beta * ones(size(eqn.U, 2), 1)];
        end
    end

    % solve the Lyapunov equation
    [adiout, eqn, opts, oper] = mess_lradi(eqn, opts, oper);

    % Store only necessary information about ADI iteration.
    tmp = adiout;

    if not(opts.nm.store_solfac)
        if isfield(tmp, 'Z'), tmp = rmfield(tmp, 'Z'); end
        if isfield(tmp, 'D'), tmp = rmfield(tmp, 'D'); end
        if isfield(tmp, 'S'), tmp = rmfield(tmp, 'S'); end
    end

    if not(opts.nm.store_debug)
        if isfield(tmp, 'DeltaK'), tmp = rmfield(tmp, 'DeltaK'); end
        if isfield(tmp, 'Knew'), tmp = rmfield(tmp, 'Knew'); end
        if isfield(tmp, 'res_fact'), tmp = rmfield(tmp, 'res_fact'); end
    end
    out.adi(k) = tmp;

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Restart Newton iteration with exact ADI iteration if necessary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if adiout.restart && opts.nm.inexact(1)

        if already_restarted
            error('MESS:lrnm', ...
                  'Newton iteration with line search failed to converge.');
        else
            % restart Newton iteration
            warning('MESS:lrnm', ...
                    ['Newton iteration needs to be restarted because ', ...
                     'of divergence. Continuing with exact ADI iteration.']);
            if opts.nm.res_tol
                res = [res(1 : k - 1); zeros(opts.nm.maxiter,1)];
            else
                res = [];
            end

            if opts.nm.rel_diff_tol
                rc = [rc(1 : k - 1); zeros(opts.nm.maxiter,1)];
            else
                rc = [];
            end

            opts.nm.maxiter = opts.nm.maxiter + k;

            % Reset U and V to initial state.
            if eqn.type == 'T'
                eqn.U = [eqn.U(:, 1:eqn.sizeUV1), -eqn.B];
                eqn.V = [eqn.V(:, 1:eqn.sizeUV1), zeros(size(eqn.B))];
            else
                eqn.U = [eqn.U(:, 1:eqn.sizeUV1), ...
                    zeros(size(eqn.C,2), size(eqn.C,1))];
                eqn.V = [eqn.V(:, 1:eqn.sizeUV1), -eqn.C'];
            end

            eqn.haveUV = 1;

            if isfield(opts.nm,'K0')
                if eqn.type == 'T'
                    eqn.V(:, end-m+1:end) = opts.nm.K0';
                else
                    eqn.U(:, end-m+1:end) = opts.nm.K0';
                end
            end

            if bdf
                if eqn.type == 'T'
                    eqn.U = (-opts.bdf.tau * opts.bdf.beta) * eqn.B;
                else
                    eqn.V = (-opts.bdf.tau * opts.bdf.beta) * eqn.C';
                end
            end

            if eqn.type == 'T'
                if eqn.diagonalized_RHS
                    eqn.G = eqn.C' * eqn.U_diag;
                else
                    eqn.G = eqn.C';
                end
            else
                if eqn.diagonalized_RHS
                    eqn.G = eqn.B * eqn.U_diag;
                else
                    eqn.G = eqn.B;
                end
            end

            eqn.G = oper.init_res(eqn, opts, oper, eqn.G);

            if opts.LDL_T
                eqn.S_diag = S;
                res0 = riccati_LR(eqn.G, [], opts, diag(eqn.S_diag), []);
                if eqn.type == 'T'
                    eqn.S_diag = [eqn.S_diag;...
                                  tau_beta * ones(size(eqn.V, 2), 1)];
                else
                    eqn.S_diag = [eqn.S_diag;...
                                  tau_beta * ones(size(eqn.U, 2), 1)];
                end
            else
                res0 = norm(eqn.G'*eqn.G, opts.norm);
                eqn.S_diag = [];
            end

            opts.nm.res0     = res0;
            opts.nm.inexact  = false;
            opts.adi.inexact = false;

            if eqn.type == 'T'
                eqn.G = [eqn.G, eqn.V(:, eqn.sizeUV1+1:end)];
            else
                eqn.G = [eqn.G, eqn.U(:, eqn.sizeUV1+1:end)];
            end

            if opts.nm.linesearch
                linesearch = 0;
                W_old = eqn.G;
%               if not(isfield(eqn,'K0')) % FIXME: dead code
                    DeltaK_old = [];
%                 else
%                     %         DeltaK_old = -eqn.V( : , s + 1 : end);
%                     %         careful with haveUV option!
%                     DeltaK_old = [];
%                end
                if opts.LDL_T
                    S_old = eqn.S_diag;
                end
            end
            restarted = 1;
            continue
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform projection update if desired
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if project && not(mod(k,opts.nm.projection.freq)) && ...
       not(opts.nm.accumulateRes && ...
           (adiout.Riccati_res < opts.nm.res_tol))

        projected = 1;

        opts.nm.linesearch = 0; % no line search after projection

        if opts.LDL_T
            [adiout.Z, adiout.D, eqn.S_diag] = ...
                mess_galerkin_projection_acceleration(adiout.Z, ...
                'CARE' ,eqn, oper, opts, adiout.D);
        else
            adiout.Z = mess_galerkin_projection_acceleration(adiout.Z, ...
                'CARE', eqn, oper, opts);
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update feedback approximation if it has not been accumulated during
    % the inner iteration or after projection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if not(opts.adi.accumulateK) || projected

        if not(opts.adi.accumulateDeltaK) || projected

            [adiout, eqn, opts, oper ] = ...
                mess_accumulateK(eqn, opts, oper, adiout, [], adiout.Z);
        else
            if eqn.type == 'T'
                adiout.Knew = eqn.V(:, end-m+1:end) + adiout.DeltaK;
            else
                adiout.Knew = eqn.U(:, end-m+1:end) + adiout.DeltaK;
            end
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.nm.res_tol||opts.nm.rel_diff_tol
        if eqn.type == 'T'
            V1 = adiout.Knew - eqn.V(:, end-m+1:end);
        else
            V1 = adiout.Knew - eqn.U(:, end-m+1:end);
        end
    end

    if opts.nm.res_tol
        if projected
            if opts.LDL_T
                res(k) = mess_res2_norms(adiout.Z,'riccati',eqn,opts,oper, ...
                         opts.nm, adiout.D)/res0;
            else
                res(k) = mess_res2_norms(adiout.Z,'riccati',eqn,opts,oper, ...
                         opts.nm,[])/res0;
            end
        else
            if opts.nm.accumulateRes
                res(k) = adiout.Riccati_res;
            else
                res(k) = riccati_LR(adiout.res_fact, V1, opts,...
                         diag(eqn.S_diag), [])/res0;
            end
        end
    end

    if opts.nm.linesearch
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check whether line search is necessary
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ((k == 1) && (res(k) > 1)) || ((k > 1) && (res(k) > res(k - 1))) ...
           || adiout.linesearch

            linesearch = 1;

            % Compute Lambda
            if opts.LDL_T
                lambda = exact_line_search(W_old, DeltaK_old, ...
                         adiout.res_fact, adiout.DeltaK, eqn.S_diag, S_old);
            else
                lambda = exact_line_search(W_old, DeltaK_old, ...
                         adiout.res_fact, adiout.DeltaK, [], []);
            end

            if opts.nm.info
                fprintf(['\n\t\t Using line search (res: %4d)\n',...
                         '\t\t lambda: %e\n'], res(k), lambda);
            end

            % Update K
            if eqn.type == 'T'
                adiout.Knew = eqn.V(:, end-m+1:end) + lambda * adiout.DeltaK;
            else
                adiout.Knew = eqn.U(:, end-m+1:end) + lambda * adiout.DeltaK;
            end

            % Update DeltaK and W
            if lambda <= 1
                adiout.DeltaK = [sqrt(1 - lambda) * DeltaK_old, ...
                                 lambda * adiout.DeltaK];
                adiout.res_fact = [sqrt(1 - lambda) * W_old, ...
                                   sqrt(lambda) * adiout.res_fact];
                if opts.LDL_T
                    S_linesearch = [S_old; eqn.S_diag];
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
                    if bdf
                        S_linesearch = [(opts.bdf.tau * opts.bdf.beta) ...
                                        * ones(size(DeltaK_old, 2), 1); S_old; ...
                                        (opts.bdf.tau * opts.bdf.beta) ...
                                        * ones(size(DeltaK_new, 2), 1)];
                        S_K = [S_old; eqn.S_diag; ...
                               (opts.bdf.tau * opts.bdf.beta) ...
                               * ones(size(DeltaK_old, 2), 1)];
                    else
                        S_linesearch = [ones(size(DeltaK_old, 2), 1); S_old; ...
                                        ones(size(DeltaK_new, 2), 1)];
                        S_K = [S_old; eqn.S_diag; ones(size(DeltaK_old, 2), 1)];
                    end
                end
            end

            if not(opts.LDL_T)
                S_linesearch = eqn.S_diag;
                S_K = [];
            end

            % Compute residual norm after line search
            res(k) = riccati_LR(adiout.res_fact, adiout.DeltaK, opts, ...
                     diag(S_linesearch), diag(S_K)) / res0;

            if not(restarted) && ...
               (((k == 1) && (res(k) >= (1 - lambda * alpha))) || ...
                ((k > 1) && (res(k) >= (1 - lambda * alpha) * res(k - 1))))

                % No sufficient decrease
                warning('MESS:lrnm', ...
                        ['Newton iteration with line search has', ...
                         ' insufficient decrease. ']);
                if opts.adi.inexact
                    % switch to exact ADI iteration
                    warning('MESS:lrnm', ...
                            'Switching to exact ADI iteration.');
                    opts.adi.inexact = 0;
                else
                    % Newton iteration with line search failed to converge, stop
                    error('MESS:lrnm', ...
                          ['Newton iteration with line search failed ' ...
                           'to converge.']);
                end
            end
        end

        % Keep ADI LR residual factor and DeltaK for next Newton iteration
        W_old = adiout.res_fact;
        DeltaK_old = adiout.DeltaK;

        if opts.LDL_T
            if linesearch
                S_old = S_linesearch;
                % be careful here, adiout.res_fact, adiout.DeltaK don't fit
                % out.S anymore, till now they are not needed anywhere
                % together after this point
            else
                S_old = eqn.S_diag;
            end
        end

        linesearch = 0;
    end

    if opts.nm.rel_diff_tol
        if opts.adi.accumulateDeltaK
            if projected
                if eqn.type == 'T'
                    adiout.DeltaK = adiout.Knew - eqn.V(:, end-m+1:end);
                else
                    adiout.DeltaK = adiout.Knew - eqn.U(:, end-m+1:end);
                end
            end
            rc(k) = norm(adiout.DeltaK, opts.norm) ...
                    / norm(adiout.Knew, opts.norm);
        else
            rc(k) = norm(V1, opts.norm);
        end
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.nm.info
        if opts.nm.rel_diff_tol&&opts.nm.res_tol
            fprintf(1,['\n NM step: %4d  normalized residual: %e\n' ...
                       '               relative change in K: %e\n' ...
                       '               number of ADI steps: %d \n\n'], ...
                    k,res(k),rc(k),adiout.niter);
        elseif opts.nm.res_tol
            fprintf(1,['\n NM step: %4d  normalized residual: %e \n\n' ...
                       '                number of ADI steps: %d \n\n'], ...
                    k,res(k),adiout.niter);
        elseif opts.nm.rel_diff_tol
            fprintf(1,['\n NM step: %4d  relative change in K: %e\n\n' ...
                       '               number of ADI steps: %d \n\n'], ...
                    k,rc(k));
        end
    end

    if eqn.type == 'T'
        eqn.V(:, end-m+1:end) = adiout.Knew;
    else
        eqn.U(:, end-m+1:end) = adiout.Knew;
    end

    k = k + 1;

    if restarted
        already_restarted = 1;
        restarted = 0;
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (opts.nm.res_tol && (res(k-1) < opts.nm.res_tol)) || ...
       (opts.nm.rel_diff_tol && (rc(k-1) < opts.nm.rel_diff_tol))
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
                opts.adi.inexact = res(k - 1) < 1 / (k ^ 3 + 1) * res(k - 2);
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
            warning('MESS:lrnm','Turning inexact ADI iteration back on.');
        end
    end

    if opts.adi.inexact
        switch opts.nm.inexact
            case 'linear'
                opts.adi.outer_tol = opts.nm.tau * res(k - 1);
            case 'superlinear'
                opts.adi.outer_tol = 1 / (k ^ 3 + 1) * res(k - 1);
            case 'quadratic'
                if res(k - 1) > 1
                    opts.adi.outer_tol = opts.nm.tau / sqrt(res(k - 1));
                else
                    opts.adi.outer_tol = opts.nm.tau * res(k - 1) * res(k - 1);
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
    out.K = eqn.V(:, end-m+1:end)';
else
    out.K = eqn.U(:, end-m+1:end)';
end
if opts.LDL_T && isfield(adiout, 'Z') && not(isempty(adiout.Z))
    out.D = adiout.D;
    out.Z = adiout.Z;
else
    if isfield(adiout, 'Z') && not(isempty(adiout.Z))
        out.Z = mess_column_compression(adiout.Z, 'N', [], ...
                opts.nm.trunc_tol, opts.nm.trunc_info);
    end
end

eqn = rmfield(eqn, 'S_diag');

if opts.nm.res_tol, out.res = res(1:out.niter); end
if opts.nm.rel_diff_tol, out.rc = rc(1:out.niter); end

% Delete the res0 part from opts for later use of the struct.
if isfield(opts.nm, 'res0')
    out.res0 = opts.nm.res0;
    opts.nm  = rmfield(opts.nm, 'res0');
elseif exist('res0', 'var')
    out.res0 = res0;
end

if out.niter == opts.nm.maxiter
    warning('MESS:NM:convergence', ...
            ['LR-NM reached maximum iteration number.', ...
             'Results may be inaccurate!']);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (size(eqn.V, 2) > eqn.sizeUV1) || (size(eqn.U, 2) > eqn.sizeUV1)
    % Cut off the stabilizing feedback.
    eqn.V = eqn.V(:, 1:eqn.sizeUV1);
    eqn.U = eqn.U(:, 1:eqn.sizeUV1);
end

if isempty(eqn.V) || isempty(eqn.U)
    % Enforce empty matrices and parameters.
    eqn.U       = [];
    eqn.V       = [];
    eqn.haveUV  = 0;
    eqn.sizeUV1 = 0;
end

% Delete overwritten right hand-side.
eqn = rmfield(eqn, 'G');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[eqn,opts,oper] = oper.mul_E_post(eqn,opts,oper);

[eqn, opts, oper] = oper.init_res_post(eqn, opts, oper);
