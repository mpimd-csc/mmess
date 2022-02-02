function [out, eqn, opts, oper] = mess_lradi(eqn, opts, oper)

%% function [out, eqn, opts, oper] = mess_lradi(eqn, opts, oper)
%
% Solve continuous-time Lyapunov equations with sparse coefficients
%   eqn.type = 'N' -> A*Z*Z'*E' + E*Z*Z'*A' + B*B' = 0 (N)
%   eqn.type = 'T' -> A'*Z*Z'*E + E'*Z*Z'*A + C'*C = 0 (T)
%
%
% Matrix A can have the form A = Ã + U*V' if U (eqn.U) and V (eqn.V) are
% provided U and V are dense (n x m3) matrices and should satisfy m3 << n
%
% Input
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation
%                       with A and E
%
% Output
%   out                 struct containing output information
%
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation
%                       with A and E
%
% Input fields in struct eqn:
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   eqn.G       dense (n x m1) matrix G
%               if present it is used instead of B, or C' as RHS
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
%               (optional, default fallback: 'N')
%
%   eqn.haveE   possible  values: 0, 1, false, true
%               if haveE = 0: matrix E is assumed to be the identity
%               (optional, default: 0)
%
%   eqn.haveUV  possible  values: 0, 1, false, true
%               if haveUV = 1: U = [U1, U2] and V = [V1, V2]
%               if K or DeltaK are accumulated during the iteration they
%               use only U2 and V2. U1 and V1 can be used for an external
%               rank-k update of the operator.
%               The size of U1 and V1 can be given via eqn.sizeUV1.
%               (optional, default: 0 if no U and V are given)
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
%   opts.norm                   possible  values: 2, 'fro'
%                               use 2-norm (2) or Frobenius norm ('fro') to
%                               compute residual and relative change norms
%                               (optional, default: 'fro')
%
%   opts.LDL_T                  possible  values: 0, 1, false, true
%                               use LDL^T formulation for the RHS and
%                               solution
%                               (optional, default: 0)
%
%   opts.adi.maxiter            possible  values: integer > 0
%                               maximum iteration number
%                               (optional, default: 100)
%
%   opts.adi.res_tol             possible  values: scalar >= 0
%                               stopping tolerance for the relative
%                               residual norm; if res_tol = 0 the relative
%                               residual norm is not evaluated
%                               (optional, default: 0)
%
%   opts.adi.rel_diff_tol              possible  values: scalar >= 0
%                               stopping tolerance for the relative
%                               change of the solution Z;
%                               if res_tol = 0 the relative
%                               change is not evaluated
%                               (optional, default: 0)
%
%   opts.adi.info               possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the status output in
%                               every iteration step
%                               (optional, default: 0)
%
%   opts.adi.compute_sol_fac    possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the computation of
%                               the factored solution; turn off if only the
%                               feedback matrix K is of interest
%                               (optional, default: 1)
%
%   opts.adi.accumulateK        possible  values: 0, 1, false, true
%                               accumulate the feedback matrix K during the
%                               iteration
%                               (optional, default: 0)
%
%   opts.adi.accumulateDeltaK   possible  values: 0, 1, false, true
%                               accumulate the update DeltaK of the
%                               feedback matrix K during the iteration
%                               (optional, default: 0)
%
%   opts.shifts.p               array with ADI shifts
%                               complex shifts are possible
%                               (optional if opts.shifts.method =
%                               'projection')
%
%   opts.shifts.info            possible  values: 0, 1, false, true
%                               turn output of used shifts before the first
%                               iteration step on (1) or off (0)
%                               (optional, default: 0)
%
% For the following fields see mess_para:
%
%   opts.shifts.method
%
%
% If optional input arguments are missing they may be set to default values
% and a often a 'MESS:control_data' warning is printed. To turn these
% warnings off use warning('OFF', 'MESS:control_data')
%
% Matrix A can have the form A = Ã + U*V'
%     if U (eqn.U) and V (eqn.V) are provided
%     U and V are dense (n x m3) matrices and should have full rank m3 << n
%     in solvers they will be treated by the Sherman-Morrison-Woodbury
%     formula.
%
% When used as the inner method, e.g. in a Newton-Kleinman method,
% the feedback matrix K can be accumulated during the iteration:
%     eqn.type = 'N' -> K = C ZZ' E
%     eqn.type = 'T' -> K = B' ZZ' E
%
% For LDL^T formulation use opts.LDL_T = 1:
%     A*L*D*L'*E' + E*L*D*L'*A' + G*S*G' = 0
%     RHS has form G * S * G'
%     Solution has form L * D * L'
%     L is stored in Z if computed (opts.adi.compute_sol_fac)
%     G (eqn.G) and S (eqn.S) need to be given
%
% Output fields in struct out:
%   out.Z               low rank solution factor
%
%   out.S               vector with diagonal elements of diagonalized
%                       eqn.S = U * out.S * U'; U is multiplied to the RHS
%
%   out.D               solution factor for LDL^T formulation
%                       (opts.LDL_T = 1)
%
%   out.res             array of relative residual norms
%
%   out.rc              array of relative change norms
%
%   out.niter           number of ADI iterations
%
%   out.res_fact        low rank residual factor W
%
%   out.Riccati_res     outer Riccati residual norm for Newton iteration
%                       (opts.nm.accumulateRes = 1)
%
%   out.linesearch      flag to trigger line search in Newton iteration
%                       (opts.adi.inexact ~= 0)
%
%   out.restart         flag to trigger complete restart of Newton
%                       iteration because of divergence
%
% uses operator functions size, init, init_res, init_res_pre, init_res_post,
% init_res_post, sol_ApE, mul_E, mul_E_pre, mul_E_post
%
%   See also mess_para, operatormanager.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check field opts.adi
if not(isfield(opts, 'adi')) || not(isstruct(opts.adi))
    error('MESS:control_data', ['No adi control data found in options', ...
          'structure.']);
end

%% check field opts.shifts
if not(isfield(opts, 'shifts')) || not(isstruct(opts.adi))
    warning('MESS:control_data', ['No shift computation control data ', ...
            'found in options structure.']);
    opts.shifts.info = 0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check info parameter for output verbosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.adi, 'info'))
    opts.adi.info = 0;
else
    if not(isnumeric(opts.adi.info)) && not(islogical(opts.adi.info))
        error('MESS:info', ...
              'opts.adi.info parameter must be logical or numeric.');
    end
end

if not(isfield(opts.shifts, 'info'))
    opts.shifts.info = 0;
else
    if not(isnumeric(opts.shifts.info)) && not(islogical(opts.shifts.info))
        error('MESS:info', ...
              'opts.shifts.info parameter must be logical or numeric.');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check stopping parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.adi, 'maxiter')) || not(isnumeric(opts.adi.maxiter))
    warning('MESS:control_data', ...
            ['Missing or Corrupted opts.adi.maxiter field. ', ...
             'Switching to default: 100']);
    opts.adi.maxiter = 100;
end

if not(isfield(opts.adi, 'rel_diff_tol')) || ...
   not(isnumeric(opts.adi.rel_diff_tol))
    warning('MESS:control_data', ...
            ['Missing or Corrupted opts.adi.rel_diff_tol field. ', ...
             'Switching to default: 0']);
    opts.adi.rel_diff_tol = 0;
end

if opts.adi.rel_diff_tol
    nrmZ = 0;
end

if not(isfield(opts.adi, 'res_tol')) || not(isnumeric(opts.adi.res_tol))
    warning('MESS:control_data', ...
            ['Missing or Corrupted opts.adi.res_tol field. ', ...
             'Switching to default: 0']);
    opts.adi.res_tol = 0;
end

if not(isfield(opts, 'norm')) || ...
   (not(strcmp(opts.norm, 'fro')) && ...
    (not(isnumeric(opts.norm)) || opts.norm ~= 2))

    warning('MESS:control_data', ...
            ['Missing or Corrupted opts.norm field. ', ...
             'Switching to default: ''fro''']);
    opts.norm = 'fro';
end

if not(isfield(opts.adi, 'inexact')), opts.adi.inexact = 0; end

if opts.adi.inexact
    if not(opts.adi.res_tol)
        % res_tol is needed
        opts.adi.res_tol = 1e-16;
        opts.adi.accumulateDeltaK = 1;
    end

    if not(isfield(opts.adi, 'outer_tol'))
        error('MESS:outer_tol', ...
              'For inexact ADI opts.adi.outer_tol is needed.');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'type'))
    eqn.type = 'N';
    warning('MESS:control_data', ['Unable to determine type of ', ...
                                  'equation. Falling back to type ''N''']);
elseif (eqn.type ~= 'N') && (eqn.type ~= 'T')
    error('MESS:equation_type', ...
          'Equation type must be either ''T'' or ''N''');
end

%set flag 0 if E does not exist
if not(isfield(eqn, 'haveE'))
    eqn.haveE = 0;
    warning('MESS:control_data', ...
            ['Missing or Corrupted eqn.haveE field.', ...
             'Switching to default: 0']);
end

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');

if not(result)
    error('MESS:control_data', ...
          'system data is not completely defined or corrupted');
end

if eqn.type == 'N' && ...
   (isfield(opts.adi, 'accumulateDeltaK') && opts.adi.accumulateDeltaK)
    if not(isfield(eqn, 'B')) || not(isnumeric(eqn.B))
        error('MESS:control_data', 'eqn.B is not defined or corrupted');
    end

    if not(isfield(eqn, 'C')) || not(isnumeric(eqn.C))
        error('MESS:control_data', 'eqn.C is not defined or corrupted');
    end

    m = size(eqn.C, 1);
end

if eqn.type == 'T' && ...
   (isfield(opts.adi, 'accumulateDeltaK') && opts.adi.accumulateDeltaK)
    if not(isfield(eqn, 'C')) || not(isnumeric(eqn.C))
        error('MESS:control_data', 'eqn.C is not defined or corrupted');
    end

    if not(isfield(eqn, 'B')) || not(isnumeric(eqn.B))
        error('MESS:control_data', 'eqn.B is not defined or corrupted');
    end

    m = size(eqn.B, 2);
end

% make sure the first right hand side is dense so that the resulting factor
% is densely stored.
if isfield(eqn, 'G') && issparse(eqn.G), eqn.G = full(eqn.G); end
if isfield(eqn, 'B') && issparse(eqn.B), eqn.B = full(eqn.B); end
if isfield(eqn, 'C') && issparse(eqn.C), eqn.C = full(eqn.C); end
if isfield(eqn, 'U') && issparse(eqn.U), eqn.U = full(eqn.U); end
if isfield(eqn, 'V') && issparse(eqn.V), eqn.V = full(eqn.V); end

% check whether LDL^T formulation should be used
if not(isfield(opts, 'LDL_T')), opts.LDL_T = 0; end

% check for or set proper right hand side in eqn.G
if opts.LDL_T
    % RHS has form G * S * G'
    % Solution has form L * D * L' with D Kronecker product of out.D and S
    % D is not computed explicitly
    % L is stored in Z if computed (opts.adi.compute_sol_fac)
    % G (eqn.G) and S (eqn.S) need to be given
    if not(isfield(eqn, 'G')) || not(isnumeric(eqn.G))
        error('MESS:control_data', 'eqn.G is not defined or corrupted');
    end

    if not(isfield(eqn, 'S')) || not(isnumeric(eqn.S))
        error('MESS:control_data', 'eqn.S is not defined or corrupted');
    end

    % init solution factor D
    out.D = zeros(opts.adi.maxiter, opts.adi.maxiter);

    if isfield(eqn, 'S_diag')
        diagonalized_RHS = 0;
    elseif isdiag(eqn.S) %%% enq.S can be a vector from lrnm and then
        % U is unknown
        eqn.S_diag = diag(eqn.S);
        diagonalized_RHS = 0;
    else
        % diagonalze S and use U to transform the initial RHS later
        [eqn.U_diag, eqn.S_diag] = eig(eqn.S);
        eqn.S_diag = diag(eqn.S_diag);
        diagonalized_RHS = 1;
        eqn.diagonalized_RHS = 1;
    end
else
    diagonalized_RHS = 0;
    if not(isfield(eqn, 'G'))
        if eqn.type == 'N'
            eqn.G = eqn.B;
        else
            eqn.G = eqn.C';
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for shifts and their properness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
init_shifts = 0;

if not(isfield(opts, 'shifts')) || not(isstruct(opts.shifts))
    error('MESS:control_data', ...
          'shift parameter control structure missing.');
end

if isfield(opts.shifts, 'method') && ...
   strcmp(opts.shifts.method, 'projection')

    opts.adi.compute_sol_fac = 1;
    opts.shifts.used_shifts = [];

    if not(isfield(opts.shifts, 'p'))
        init_shifts = 1;
    end

    if not(isfield(opts.shifts, 'num_desired'))
        if opts.LDL_T
            opts.shifts.num_desired = max(5, size(eqn.G, 2));
        elseif eqn.type == 'N'
            opts.shifts.num_desired = max(5, size(eqn.B, 2));
        else
            opts.shifts.num_desired = max(5, size(eqn.C, 1));
        end
    end
else
    if not(isfield(opts.shifts, 'p'))
        init_shifts = 1;
    else
        illegal_shifts = 0;
        % Check if all shifts are in the open left half plane
        if any(not((real(opts.shifts.p)) < 0)), illegal_shifts = 1; end

        % Check if complex pairs of shifts are properly ordered.
        k = 1;
        while k <= length(opts.shifts.p)
            if not((isreal(opts.shifts.p(k))))
                if (opts.shifts.p(k+1) ~= conj(opts.shifts.p(k)))
                    illegal_shifts = 1;
                end
                k = k + 1;
            end
            k = k + 1;
        end
        if illegal_shifts
            error('MESS:shifts_improper', ...
                  'Improper shift vector detected!');
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for feedback and shift matrices appearing inside
% Newton, BDF and Rosenbrock type methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts, 'rosenbrock')), opts.rosenbrock = []; end

if isstruct(opts.rosenbrock) && isfield(opts.rosenbrock, 'tau')
    rosenbrock = 1;
else
    rosenbrock = 0;
end

if not(isfield(opts, 'bdf')), opts.bdf = []; end

if isstruct(opts.bdf) && isfield(opts.bdf, 'tau') && isfield(opts.bdf, 'beta')
    bdf = 1;
else
    bdf = 0;
end

% Check for rank-k update of the operator.
if not(isfield(eqn, 'U')) || isempty(eqn.U) || ...
   not(isfield(eqn, 'V')) || isempty(eqn.V)
    eqn.haveUV = 0;
else
    if isnumeric(eqn.U) && isnumeric(eqn.V) && ...
       size(eqn.U, 1) == size(eqn.V, 1) && size(eqn.U, 2) == size(eqn.V, 2)
        eqn.haveUV = 1;
    else
        error('MESS:control_data', ...
              ['Inappropriate data of low rank updated operator ', ...
               '(eqn.U and eqn.V)']);
    end
end

% Check for size of constant term in U and V.
if eqn.haveUV
    if not(isfield(eqn, 'sizeUV1')) || isempty(eqn.sizeUV1)
        eqn.sizeUV1 = size(eqn.U, 2);
    else
        assert(isnumeric(eqn.sizeUV1) && (eqn.sizeUV1 <= size(eqn.U, 2)), ...
               'MESS:control_data', ...
               ['Inappropriate size of low rank updated operator ', ...
                '(eqn.U and eqn.V)']);
    end
else
    eqn.sizeUV1 = 0;
end

% Get sizes of right hand side
k = size(eqn.G, 2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether we want to compute Z or rather accumulate K in the Newton
% context for AREs, or both, e.g., in inexact Newton contexts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accumulation of K or DeltaK is helpful in inexact Newton and implicit
% Newton settings
if isfield(opts.adi, 'accumulateK') && opts.adi.accumulateK
    if eqn.type == 'T'
        out.Knew = zeros(size(eqn.B));
    else
        out.Knew = zeros([size(eqn.C, 2), size(eqn.C, 1)]);
    end
else
    opts.adi.accumulateK = 0;
end

if isfield(opts.adi, 'accumulateDeltaK') && opts.adi.accumulateDeltaK
    if eqn.type == 'T'
        if eqn.haveUV && not(eqn.sizeUV1)
            % eqn.V is only K
            out.DeltaK = -eqn.V;
        elseif eqn.haveUV && (size(eqn.V, 2) > m)
            % eqn.V is given and K is in the second part.
            out.DeltaK = -eqn.V(:, end-m+1:end);
        else
            % K = []
            out.DeltaK = zeros(size(eqn.B));
        end
    else
        if eqn.haveUV && not(eqn.sizeUV1)
            % eqn.U is only K
            out.DeltaK = -eqn.U;
        elseif eqn.haveUV && (size(eqn.U, 2) > m)
            % eqn.U is given and K is in the second part.
            out.DeltaK = -eqn.U(:, end-m+1:end);
        else
            % K = []
            out.DeltaK = zeros([size(eqn.C, 2), size(eqn.C, 1)]);
        end
    end
else
    opts.adi.accumulateDeltaK = 0;
end

if not(isfield(opts.adi, 'compute_sol_fac'))
    opts.adi.compute_sol_fac = 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All checks done. Here comes the real work!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize required usf for multiplication with E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eqn.haveE, [eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper); end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.adi.compute_sol_fac
    Z = zeros(size(eqn.G, 1), opts.adi.maxiter*k);
else
    Z = [];
end

if opts.adi.res_tol
    res = zeros(1, opts.adi.maxiter);
else
    res = [];
end

if opts.adi.rel_diff_tol
    rc = zeros(1, opts.adi.maxiter);
else
    rc = [];
end

[eqn, opts, oper] = oper.init_res_pre(eqn, opts, oper);
% in the LDL_T case we may have diagonalized the kernel matrix of the RHS.
% If so, we need to initialize the residual with the updated G matrix
if diagonalized_RHS
    [W, res0, eqn, opts, oper] = ...
        oper.init_res(eqn, opts, oper, eqn.G*eqn.U_diag);
else
    [W, res0, eqn, opts, oper] = ...
        oper.init_res(eqn, opts, oper, eqn.G);
end

% Initialize shift vector in case of projection shifts and empty initial
% shift vector
if init_shifts
    [p, ~, eqn, opts, oper] = mess_para(eqn, opts, oper);
    opts.shifts.p = p; % opts overrides opts.shifts.p above
end

% Get length of shift vector
l = length(opts.shifts.p);

if opts.shifts.info
    fprintf('ADI Shifts:\n');
    disp(opts.shifts.p);
end

out.linesearch = 0;
out.restart = 0;

if isfield(opts, 'nm') && isfield(opts.nm, 'accumulateRes') && ...
   opts.nm.accumulateRes && isfield(opts.nm, 'res0')

    outer_res = zeros(1, opts.adi.maxiter);
    res0 = opts.nm.res0;
else
    outer_res = [];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 1;
m_shift = 1;

while m < opts.adi.maxiter + 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check whether shifts need to be updated
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if m_shift > l
        m_shift = 1;
        if strcmp(opts.shifts.method, 'projection')
            if opts.LDL_T
                [opts, l] = mess_get_projection_shifts(eqn, opts, oper, ...
                                Z(:, 1:(m - 1)*k), W, out.D(1:m-1, 1:m-1));
            else
                [opts, l] = mess_get_projection_shifts(eqn, opts, oper, ...
                                Z(:, 1:(m - 1)*k), W);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get current shift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pc = opts.shifts.p(m_shift);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if bdf
        [V, eqn, opts, oper] = ...
            mess_solve_shifted_system_BDF(eqn, opts, oper, pc, W);
    elseif rosenbrock
        [V, eqn, opts, oper] = ...
            mess_solve_shifted_system_Rosenbrock(eqn, opts, oper, pc, W);
    else
        [V, eqn, opts, oper] = ...
            mess_solve_shifted_system(eqn, opts, oper, pc, W);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update low rank solution factor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isreal(pc)
        % just update the factor
        V = real(V);
        if opts.adi.compute_sol_fac
            if opts.LDL_T
                Z(:, (m - 1)*k+1:m*k) = V;
                out.D(m, m) = -2.0 * pc;
            else
                Z(:, (m - 1)*k+1:m*k) = sqrt(-2.0*pc) * V;
            end
        end

        % update low rank residual
        if eqn.haveE
            EV = oper.mul_E(eqn, opts, eqn.type, V, 'N');
            W = W - 2.0 * pc * EV;
        else
            W = W - 2.0 * pc * V;
        end
        [out, eqn, opts, oper] = ...
            mess_accumulateK(eqn, opts, oper, out, pc, V);
    else
        % perform a double step with the known solution for the conjugate
        % shift
        a = 2.0 * sqrt(-real(pc));
        b = real(pc) / imag(pc);
        V1 = a * (real(V) + b * imag(V));
        V2 = (a * sqrt(b*b+1)) * imag(V);

        if opts.adi.compute_sol_fac
            if opts.LDL_T
                Z(:, (m - 1)*k+1:(m + 1)*k) = ...
                    [(sqrt(2.0) / a) * V1, (sqrt(2.0) / a) * V2];
                out.D(m : m + 1, m : m + 1) = -2.0 * real(pc) * eye(2);
            else
                Z(:, (m - 1)*k+1:(m + 1)*k) = [V1, V2];
            end
        end

        [out, eqn, opts, oper] = ...
             mess_accumulateK(eqn, opts, oper, out, pc, V1, V2);

        % update low rank residual for double step
        if eqn.haveE
            EV = oper.mul_E(eqn, opts, eqn.type, V1, 'N');
            W = W + a * EV;
        else
            W = W + a * V1;
        end

        m = m + 1;
        m_shift = m_shift + 1;

        if not(isempty(outer_res))
            if m > 2
                outer_res(m - 1) = outer_res(m - 2);
            else
                outer_res(m - 1) = opts.nm.res0;
            end
        end

        if not(isempty(res))
            if m > 2
                res(m - 1) = res(m - 2);
            else
                res(m - 1) = res0;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.adi.res_tol
        if opts.LDL_T
            if opts.norm == 2
                res(m) = max(abs(eig(W'*W*diag(eqn.S_diag)))) / res0;
            elseif strcmp(opts.norm, 'fro')
                res(m) = norm(eig(W'*W*diag(eqn.S_diag)), 'fro') / res0;
            end
        else
            res(m) = norm(W'*W, opts.norm) / res0;
        end

        if not(isempty(outer_res)) %riccati_LR does the LDL_T check itself.
            outer_res(m) = riccati_LR(W, out.DeltaK, opts, ...
                                      diag(eqn.S_diag), []) / opts.nm.res0;
        end
    end

    if opts.adi.rel_diff_tol
        if isreal(pc)
            nrmV = -2.0 * pc * sum(sum(V.^2));
        else % complex double step means 2 blocks added
            nrmV = sum(sum([V1, V2].^2));
        end

        nrmZ = nrmZ + nrmV;
        rc(m) = sqrt(nrmV/nrmZ);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.adi.info
        if opts.adi.rel_diff_tol && opts.adi.res_tol
            fprintf(1, ['ADI step: %4d normalized residual: %e ', ...
                        'relative change in Z: %e\n'], m, res(m), rc(m));
        elseif opts.adi.res_tol
            fprintf(1, 'ADI step: %4d normalized residual: %e \n', ...
                    m, res(m));
        elseif opts.adi.rel_diff_tol
            fprintf(1, ['ADI step: %4d relative change ', ...
                    'in Z: %e\n'], m, rc(m));
        end

        if not(isempty(outer_res))
            fprintf(1, '\t\t normalized outer residual: %e\n', ...
                    outer_res(m));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [opts, out, stop] = prepare_next_adi_iteration(opts, out, res, ...
        rc, outer_res, m);

    if stop
        break
    end

    m = m + 1;
    m_shift = m_shift + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print outer tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.adi.info && opts.adi.inexact
    fprintf(1, '\n outer tolerance: %e\n', opts.adi.outer_tol);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare output arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.niter = m - (m > opts.adi.maxiter);

if opts.adi.compute_sol_fac
    out.Z = Z(:, 1:out.niter*k);
    if opts.LDL_T
        out.D = kron(out.D(1:out.niter, 1:out.niter), diag(eqn.S_diag));
    end
end

if opts.adi.res_tol, out.res = res(1:out.niter); end

if opts.adi.rel_diff_tol, out.rc = rc(1:out.niter); end

out.res_fact = W;

if opts.LDL_T
    out.S = eqn.S_diag;
end

if not(isempty(outer_res))
    out.Riccati_res = outer_res(out.niter);
end

if out.niter == opts.adi.maxiter
    warning('MESS:ADI:convergence',...
            ['LR-ADI reached maximum iteration number.',...
             'results may be inaccurate!']);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalize required usf for multiplication with E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eqn.haveE, [eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper); end

[eqn, opts, oper] = oper.init_res_post(eqn, opts, oper);
