function [out, eqn, opts, oper] = mess_lrradi(eqn, opts, oper)
%% function [out, eqn, opts, oper] = mess_lrradi(eqn,opts, oper)
%
% Solve continuous-time Riccati equations with sparse coefficients with
% the RADI method [1]. With X = Z*inv(Y)*Z',
%   eqn.type = 'N' -> A*X*E' + E*X*A' - E*X*C'*C*X*E' + B*B' = 0 (N)
%   eqn.type = 'T' -> A'*X*E + E'*X*A - E'*X*B*B'*X*E + C'*C = 0 (T)
%
% Matrix A can have the form A = Ã + U*V' if U (eqn.U) and V (eqn.V) are
% provided U and V are dense (n x m3) matrices and should satisfy m3 << n
%
% Input/Output
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation
%                       with A and E
%
% Output
%   out                 struct containing solutions and output information
%
% Input fields in struct eqn:
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   eqn.S       dense (m1 x m1) matrix (N) or (m2 x m2) matrix (T)
%               expected to be symmetric
%               (required for LDL^T formulation)
%
%   eqn.U       dense (n x m3) matrix U
%               (optional, required if eqn.V is present)
%
%   eqn.V       dense (n x m3) matrix V
%               (optional, required if eqn.U is present)
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional, default 'N')
%
%   eqn.haveE   possible  values: 0, 1, false, true
%               if haveE = 0: matrix E in eqn.E_ is assumed to be identity
%               (optional, default 0)
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
%   opts.radi.Z0                possible values: dense (n x m4) matrix
%                               initial stabilizing solution factor
%                               X0 = Z0*inv(Y0)*Z0', this factor has to
%                               result in a positive semi-definite Riccati
%                               residual W0
%                               (optional, default: zeros(n, m4))
%
%   opts.radi.Y0                possible values: dense (m4 x m4) matrix
%                               initial stabilizing solution factor
%                               X0 = Z0*inv(Y0)*Z0', this factor has to
%                               result in a positive semi-definite
%                               Riccati residual W0
%                               (optional, default: eye(m4))
%
%   opts.radi.W0                possible values: dense (n x m5) matrix
%                               initial Riccati residual factor such that
%                               R(X0) = W0 * W0', if 
%                               opts.radi.compute_res = 1, this factor is
%                               computed out of Z0 and Y0
%                               Note: In case of Bernoulli stabilization
%                               the W0 is given by the right hand-side C'
%                               for 'T' and B for 'N' and is automatically
%                               set if opts.radi.compute_res = 0
%                               (optional, default: C' for 'T' or B for 
%                               'N')
%
%   opts.radi.S0                possible values: dense (m5 x m5) matrix
%                               initial Riccati residual factor such that
%                               R(X0) = W0 * S0 * W0', if 
%                               opts.radi.compute_res = 1, this factor is
%                               computed out of Z0 and Y0
%                               (required for LDL^T formulation if
%                                opts.radi.W0 was explicitly set)
%
%   opts.radi.K0                possible values: dense 'T': (m1 x n) 
%                               matrix, 'N':  (m2 x n) matrix
%                               initial K (corresponding to Z0 and Y0)
%                               Note: If K0 is given without Z0, only the
%                               resulting stabilizing feedback is computed.
%                               Also it has to correspond to W0.
%                               (optional, default: E*Z0*inv(Y0)*Z0'*C' for
%                               'N' or E'*Z0*inv(Y0)*Z0'*B for 'T')
%
%   opts.radi.compute_sol_fac   possible values: 0, 1, false, true
%                               turn on (1) or off (0) to compute the 
%                               solution of the Riccati equation and use it
%                               internally for computations, or only
%                               the stabilizing feedback
%                               (optional, default: 1)
%
%   opts.radi.get_ZZt           possible values: 0, 1, false, true
%                               turn on (1) or off (0) to compute only
%                               the low-rank decomposition X = Z*Z'
%                               without the middle term Y
%                               (optional, default: 1)
%
%   opts.radi.compute_res       possible values: 0, 1, false, true
%                               turn on (1) or off (0) to compute the
%                               residual corresponding to the initial
%                               solution factors Z0, Y0, if 0 then the
%                               right hand-side is used as residual if
%                               there is no W0
%                               (optional, default: 1)
%
%   opts.radi.maxiter           possible  values: integer > 0
%                               maximum RADI iteration number
%                               (optional, default: 100)
%
%   opts.radi.res_tol           possible  values: scalar >= 0
%                               stopping tolerance for the relative
%                               RADI residual norm; if res_tol = 0 the
%                               relative residual norm is not evaluated
%                               (optional, default: 0)
%
%   opts.radi.rel_diff_tol      possible  values: scalar >= 0
%                               stopping tolerance for the relative
%                               change of the RADI solution Z;
%                               if res_tol = 0 the relative
%                               change is not evaluated
%                               (optional, default: 0)
%
%   opts.norm                   possible  values: 2, 'fro'
%                               use 2-norm (2) or Frobenius norm ('fro') to
%                               compute residual and relative change norms;
%                               must be the same as opts.norm
%                               (optional, default: 'fro')
%
%   opts.radi.info              possible  values: 0, 1, false, true
%                               turn on (1) or off (0) the status output in
%                               every RADI iteration step
%                               (optional, default: 0)
%
%   opts.radi.trunc_tol         possible values: scalar > 0
%                               tolerance for rank truncation of the
%                               low-rank solutions (aka column compression)
%                               (optional, default: eps*n)
%
%   opts.radi.trunc_info        possible values: 0, 1, false, true
%                               verbose mode for column compression
%                               (optional, default: 0)
%
%   opts.shifts.method          possible  values:
%                               'precomputed',
%                               'penzl','heur', (basic MMESS routine)
%                               'projection' (basic MMESS routine)
%                               'gen-ham-opti' (special for RADI)
%                                method for shift computation
%                               (optional, default: 'gen-ham-opti')
%
%   opts.shifts.history         possible values: integer * size(W0, 2) > 0
%                               parameter for accumulating the history
%                               of shift computations
%                               (optional, default: 6 * columns of
%                               residual)
%
%   opts.shifts.info            possible  values: 0, 1, false, true
%                               turn output of used shifts before the first
%                               iteration step on (1) or off (0)
%                               (optional, default: 0)
%
%
% If optional input arguments are missing they may be set to default values
% and a 'MESS:control_data' warning is printed. To turn warnings off use
% warning('OFF', 'MESS:control_data').
%
% The feedback matrix K can be accumulated during the iteration:
%     eqn.type = 'N' -> K = (E*X*C')'
%     eqn.type = 'T' -> K = (E'*X*B)'
%
%
% Output fields in struct out:
%   out.Z           low rank solution factor, the solution is
%                   opts.radi.get_ZZt = 0: X = Z*inv(Y)*Z'
%                   opts.radi.get_ZZt = 1: X = Z*Z'
%                   (opts.radi.compute_sol_fac = 1 and not only initial K0)
%
%   out.Y           small square solution factor, the solution is
%                   opts.radi.get_ZZt = 0: X = Z*inv(Y)*Z'
%                   (opts.radi.compute_sol_fac = 1 and not only initial K0)
%
%   out.D           solution factor for LDL^T formulation, the solution is
%                   opts.LDL_T = 1: X = Z*D*Z'
%                   (opts.LDL_T = 1)
%
%   out.K           stabilizing Feedback matrix
%
%   out.timesh      time of the overall shift computation
%
%   out.p           used shifts
%
%   out.niter       number of RADI iterations
%
%   out.res         array of relative RADI residual norms
%                   (opts.radi.res_tol ~= 0)
%
%   out.rc          array of relative RADI change norms
%                   (opts.radi.rel_diff_tol ~= 0)
%
%   out.res_fact    final Riccati residual factor W of the iteration
%
%   out.res0        norm of the normalization residual term
%
%
%   uses operator functions init and mul_E, mul_E_pre, mul_E_post,
%   mul_A, mul_A_pre, mul_A_post, init_res_pre, init_res, init_res_post,
%   size directly and further indirectly
%
% References:
% [1] P. Benner, Z. Bujanović, P. Kürschner, J. Saak, RADI: A low-rank
%     ADI-type algorithm for large scale algebraic Riccati equations, 
%     Numer. Math. 138 (2) (2018) 301–330. 
%     https://doi.org/10.1007/s00211-017-0907-5.
%
%   See also mess_lrradi_get_shifts, operatormanager.

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'haveE'))
    eqn.haveE = 0;
end

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');
if not(result)
    error('MESS:control_data', ...
          'system data is not completely defined or corrupted');
end

if not(isfield(eqn, 'B')) || not(isnumeric(eqn.B))
    error('MESS:control_data', 'eqn.B is not defined or corrupted');
end

if not(isfield(eqn, 'C')) || not(isnumeric(eqn.C))
    error('MESS:control_data', 'eqn.C is not defined or corrupted');
end

% Make sure the first right hand side is dense so that the resulting factor
% is densely stored.
if issparse(eqn.C)
    eqn.C = full(eqn.C);
end

if issparse(eqn.B)
    eqn.B = full(eqn.B);
end

if not(isfield(eqn, 'type'))
    eqn.type = 'N';
    warning('MESS:control_data',['Unable to determine type of equation.'...
            'Falling back to type ''N''']);
elseif (eqn.type ~= 'N') && (eqn.type ~= 'T')
    error('MESS:equation_type', ...
          'Equation type must be either ''T'' or ''N''');
end

if not(isfield(opts, 'LDL_T')), opts.LDL_T = 0; end

if opts.LDL_T
    if not(isfield(eqn, 'S')) || not(isnumeric(eqn.S))
        error('MESS:control_data', 'eqn.S is not defined or corrupted');
    end
end

if eqn.type == 'T'
    m      = size(eqn.B, 2); %number of inputs
    eqn.BB = eqn.B;          %set quadratic term
    eqn.CC = eqn.C';         %set right hand side
else
    m      = size(eqn.C, 1); %number of outputs
    eqn.BB = eqn.C';         %set quadratic term
    eqn.CC = eqn.B;          %set right hand side
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
               'eqn.haveUV option.']);
    end
    
    if isnumeric(eqn.U) && isnumeric(eqn.V) && ...
       size(eqn.U, 1) == size(eqn.V, 1) && size(eqn.U, 2) == size(eqn.V, 2)
        
        if issparse(eqn.V), eqn.V = full(eqn.V); end
        if issparse(eqn.U), eqn.U = full(eqn.U); end
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
% Initialize required usf for multiplications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eqn.haveE
    [eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);
end

[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for RADI Control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts, 'radi')) || not(isstruct(opts.radi))
    error('MESS:control_data', ['No radi control data found in ', ...
          'options structure.']);
end

% Check computation of Riccati solution.
if not(isfield(opts.radi, 'compute_sol_fac')) || ...
   isempty(opts.radi.compute_sol_fac)
    opts.radi.compute_sol_fac = 1;
end

% Check format of output solution.
if not(isfield(opts.radi, 'get_ZZt')) || isempty(opts.radi.get_ZZt)
    opts.radi.get_ZZt = 1;
end

% Check for residual norm.
if not(isfield(opts, 'norm')) || (not(strcmp(opts.norm, 'fro')) && ...
                                  (not(isnumeric(opts.norm)) || opts.norm ~= 2))
    warning('MESS:control_data', ...
            ['Missing or Corrupted opts.norm field.', ...
             'Switching to default: ''fro''']);
    opts.norm = 'fro';
end

if not(isfield(opts.radi,'trunc_tol'))
    opts.radi.trunc_tol = eps * oper.size(eqn, opts);
end

if not(isfield(opts.radi, 'trunc_info')), opts.radi.trunc_info = 0; end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List all currently unsupported options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(opts, 'bdf') && not(isempty(opts.bdf))
    error( 'MESS:control_data', 'Options bdf not supported.');
end

if isfield(opts, 'rosenbrock') && not(isempty(opts.rosenbrock))
    error( 'MESS:control_data', 'Options rosenbrock not supported.');
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for initial values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize projected residual factor in case of DAE
% oper.init_res assumes that opts.LDL_T is defined.
[eqn, opts, oper] = oper.init_res_pre(eqn, opts, oper);

onlyK = 0;
hasZ0 = isfield(opts.radi, 'Z0') && not(isempty(opts.radi.Z0));
hasY0 = isfield(opts.radi, 'Y0') && not(isempty(opts.radi.Y0));
hasW0 = isfield(opts.radi, 'W0') && not(isempty(opts.radi.W0));
hasS0 = isfield(opts.radi, 'S0') && not(isempty(opts.radi.S0));
hasK0 = isfield(opts.radi, 'K0') && not(isempty(opts.radi.K0));

if hasZ0
    % Stabilizing initial solution.
    Z   = opts.radi.Z0;
    nZ0 = size(Z, 2);
    
    % Initial middle term.
    if hasY0
        Y = opts.radi.Y0;
    else
        Y = eye(size(Z, 2));
    end
    
    % Check for residual computation option.
    if not(isfield(opts.radi, 'compute_res')) || isempty(opts.radi.compute_res)
        warning('MESS:control_data', ...
                ['Missing or Corrupted opts.radi.compute_res field.', ...
                 'Switching to default: 1']);
        opts.radi.compute_res = 1;
    end
    
    % Initial residual.
    if hasW0
        % Case: Initial residual is given.
        W = opts.radi.W0;
        if opts.LDL_T
            assert(hasS0, ...
                   'MESS:control_data', ...
                   'Missing or corrupted opts.radi.S0 field.');
            if isdiag(opts.radi.S0)
                eqn.S_diag = diag(opts.radi.S0);
            else
                [eqn.U_diag, eqn.S_diag] = eig(opts.radi.S0);
                eqn.S_diag               = diag(eqn.S_diag);
                W                        = W * eqn.U_diag;
            end
        end

        [W, res0, eqn, opts, oper] = oper.init_res(eqn, opts, oper, W);

    elseif opts.radi.compute_res
        % Case: Initial residual has to be computed.
        AZ = oper.mul_A(eqn, opts, eqn.type, Z, 'N');

        if eqn.haveE
            EZ = oper.mul_E(eqn, opts, eqn.type, Z, 'N');
        else
            EZ = Z;
        end
        
        if eqn.sizeUV1
            if eqn.type == 'T'
                UU = eqn.V(:, 1:eqn.sizeUV1);
                VV = EZ * (Y \ (Z' * eqn.U(:, 1:eqn.sizeUV1)));
            else
                UU = eqn.U(:, 1:eqn.sizeUV1);
                VV = EZ * (Y \ (Z' * eqn.V(:, 1:eqn.sizeUV1)));
            end
            UDV = [zeros(eqn.sizeUV1), eye(eqn.sizeUV1); ...
                   eye(eqn.sizeUV1),   zeros(eqn.sizeUV1)];
        else
            UU  = [];
            VV  = [];
            UDV = [];
        end
        
        if opts.LDL_T
            D0 = blkdiag([zeros(size(Y)), Y \ eye(size(Y)); ...
                          Y \ eye(size(Y)), zeros(size(Y))], ...
                          UDV, -eye(m), eqn.S);
        else
            D0 = blkdiag([zeros(size(Y)), Y \ eye(size(Y)); ...
                          Y \ eye(size(Y)), zeros(size(Y))], ...
                          UDV, -eye(m), eye(size(eqn.CC, 2)));
        end
        
        [G, S] = mess_column_compression( ...
                  [AZ, EZ, UU, VV, EZ * (Y \ (Z' * eqn.BB)), eqn.CC], 'N', ...
                  D0, opts.radi.trunc_tol, opts.radi.trunc_info);
        
        if opts.LDL_T || not(all(diag(S) > 0))
            if not(opts.LDL_T)
                warning('MESS:control_data', ...
                        ['The initial residual is indefinite, ' ...
                         'change to LDL^T approach!']);
                opts.LDL_T = 1;
            end
            
            W          = G;
            eqn.S_diag = diag(S);
        else
            W = G * diag(sqrt(diag(S)));
        end
        
        [W, res0, eqn, opts, oper] = oper.init_res(eqn, opts, oper, W);
    else
        % Case: Initial residual is given as right hand-side (Bernoulli).
        W = eqn.CC;
        if opts.LDL_T
            if isdiag(eqn.S)
                eqn.S_diag = diag(eqn.S);
            else
                [eqn.U_diag, eqn.S_diag] = eig(eqn.S);
                eqn.S_diag               = diag(eqn.S_diag);
                W                        = W * eqn.U_diag;
            end
        end

        [W, res0, eqn, opts, oper] = oper.init_res(eqn, opts, oper, W);
    end
    
    % Initial stabilizing feedback.
    if eqn.type == 'T'
        if isfield(opts.radi, 'K0') && not(isempty(opts.radi.K0))
            eqn.V(: , end-m+1:end) = opts.radi.K0';
        else
            if eqn.haveE
                eqn.V(: , end-m+1:end) = oper.mul_E(eqn, opts, 'T', Z, 'N') ...
                                         * (Y \ (Z' * eqn.B));
            else
                eqn.V(: , end-m+1:end) = Z * (Y \ (Z' * eqn.B));
            end
        end
    else
        if isfield(opts.radi, 'K0') && not(isempty(opts.radi.K0))
            eqn.U(: , end-m+1:end) = opts.radi.K0';
        else
            if eqn.haveE
                eqn.U(: , end-m+1:end) = oper.mul_E(eqn, opts, 'N', Z, 'N') ...
                                         * (Y \ (Z' * eqn.C'));
            else
                eqn.U(: , end-m+1:end) = Z * (Y \ (Z' * eqn.C'));
            end
        end
    end
elseif hasK0
    % Stabilizing initial feedback.
    if eqn.type == 'T'
        eqn.V(: , end-m+1:end) = opts.radi.K0';
    else
        eqn.U(: , end-m+1:end) = opts.radi.K0';
    end
    
    % Assign the corresponding residual.
    if hasW0
        W = opts.radi.W0;
        
        if opts.LDL_T
            assert(hasS0, ...
                   'MESS:control_data', ...
                   'Missing or corrupted opts.radi.S0 field.');
            if isdiag(opts.radi.S0)
                eqn.S_diag = opts.radi.S0;
            else
                [eqn.U_diag, eqn.S_diag] = eig(opts.radi.S0);
                eqn.S_diag               = diag(eqn.S_diag);
                W                        = W * eqn.U_diag;
            end
        end
    else
        W = eqn.CC;
        
        if opts.LDL_T
            if isdiag(eqn.S)
                eqn.S_diag = diag(eqn.S);
            else
                [eqn.U_diag, eqn.S_diag] = eig(eqn.S);
                eqn.S_diag               = diag(eqn.S_diag);
                W                        = W * eqn.U_diag;
            end
        end
    end
    
    [W, res0, eqn, opts, oper] = oper.init_res(eqn, opts, oper, W);
    
    % Other initial values.
    Z     = zeros(oper.size(eqn, opts), 0);
    nZ0   = 0;
    Y     = [];
    onlyK = 1;
else
    % Start with zero initial solution.
    Z         = zeros(oper.size(eqn, opts), 0);
    nZ0       = 0;
    Y         = [];
    
    W = eqn.CC;
    if opts.LDL_T
        if isdiag(eqn.S)
            eqn.S_diag = diag(eqn.S);
        else
            [eqn.U_diag, eqn.S_diag] = eig(eqn.S);
            eqn.S_diag               = diag(eqn.S_diag);
            W                        = W * eqn.U_diag;
        end
    end

    [W, res0, eqn, opts, oper] = oper.init_res(eqn, opts, oper, W);
end
p = size(W, 2); % size of residual factor.


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for shift parameter structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts, 'shifts')) || not(isstruct(opts.shifts))
    error('MESS:control_data', ...
          'shift parameter control structure missing.');
end

% Default shift method settings.
if not(isfield(opts, 'shifts')) || ...
   not(isstruct(opts.shifts)) || ...
   not(isfield(opts.shifts, 'method'))

    warning('MESS:control_data',...
            ['shift parameter control structure missing.', ...
             'Switching to default: method = gen-ham-opti, history = 6.']);
    
    opts.shifts.method  = 'gen-ham-opti';
end

% Check for shift history parameter.
if not(isfield(opts.shifts, 'history')) || isempty(opts.shifts.history)
    opts.shifts.history = 6 * p;
end

% Use heuristic penzl shifts routines from MMESS.
if strcmp(opts.shifts.method, 'penzl') || ...
   strcmp(opts.shifts.method, 'heur') || ...
   strcmp(opts.shifts.method, 'projection')

    if not(isfield(opts.shifts, 'num_desired'))
        opts.shifts.num_desired = opts.shifts.history;
    end
    
    if strcmp(opts.shifts.method, 'penzl') || ...
       strcmp(opts.shifts.method, 'heur')

        if not(isfield(opts.shifts, 'num_Ritz'))
            opts.shifts.num_Ritz = opts.shifts.history + 1;
        end
        
        if not(isfield(opts.shifts,'num_hRitz'))
            opts.shifts.num_hRitz = opts.shifts.history;
        end
    end
end

% Use provided shifts.
if strcmp(opts.shifts.method, 'precomputed')
    if not((isfield(opts.shifts, 'p')) && ...
       isnumeric(opts.shifts.p) && ...
       isvector(opts.shifts.p))

        error('MESS:shifts', ...
              'Found empty shift vector. Please provide proper shifts.');
    else
        illegal_shifts = 0;
        
        % Check if all shifts are in the open left half plane
        if any(not(real(opts.shifts.p) < 0))
            illegal_shifts = 1;
        end
        
        % Check if complex pairs of shifts are properly ordered.
        k = 1;
        while k <= length(opts.shifts.p)
            if not(isreal(opts.shifts.p(k)))
                if not(opts.shifts.p(k+1) == conj(opts.shifts.p(k)))
                    illegal_shifts = 1;
                end
                k = k+1;
            end
            k = k+1;
        end
        
        if illegal_shifts
            error('MESS:shifts_improper', 'Improper shift vector detected!');
        end
    end
else
    % If the shifts are not precomputed, let the shift array initially be
    % empty.
    opts.shifts.p = [];
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check info parameter for output verbosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.radi, 'info'))
    opts.radi.info = 0;
else
    if not(isnumeric(opts.radi.info)) && not(islogical( opts.radi.info ))
        error('MESS:info', ...
              'opts.radi.info parameter must be logical or numeric.');
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
if not(isfield(opts.radi, 'maxiter')) || not(isnumeric(opts.radi.maxiter))
    warning('MESS:control_data',...
            ['Missing or Corrupted opts.radi.maxiter field.', ...
             'Switching to default: 100']);
    opts.radi.maxiter = 100;
end

if not(isfield(opts.radi,'rel_diff_tol')) || not(isnumeric(opts.radi.rel_diff_tol))
    warning('MESS:control_data',...
            ['Missing or Corrupted opts.radi.rel_diff_tol field.', ...
             'Switching to default: 0']);
    opts.radi.rel_diff_tol = 0;
end

if opts.radi.rel_diff_tol
    nrmZ = sum(sum(Z.^2));
end

if not(isfield(opts.radi, 'res_tol')) || not(isnumeric(opts.radi.res_tol))
    warning('MESS:control_data',...
            ['Missing or Corrupted opts.radi.res_tol field.', ...
             'Switching to default: 0'] );
    opts.radi.res_tol = 0;
end

% Check if the low-rank factor Z needs to be computed entirely.
maxcolZ = (opts.radi.maxiter + 1) * p + nZ0;
opts.radi.compute_sol_facpart = 0;

if (strcmp(opts.shifts.method, 'gen-ham-opti') || ...
    strcmp(opts.shifts.method, 'projection')) && ...
   (opts.radi.compute_sol_fac == 0)

    opts.radi.compute_sol_facpart = 1;
    maxcolZ = opts.shifts.history;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All checks done. Here comes the real work!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.radi.res_tol
    res = zeros(1, opts.radi.maxiter);
else
    res = [];
end

if opts.radi.rel_diff_tol
    rc = zeros(1, opts.radi.maxiter);
else
    rc = [];
end

% Get relevant sizes of right hand side and shift vector.
nShifts = length(opts.shifts.p);    

if opts.shifts.info  % print shifts
    fprintf('RADI Shifts:\n');
    disp(opts.shifts.p);
end

% Reset the timer.
out.timesh = 0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
k_shift = 1;

while k < (opts.radi.maxiter + 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check whether shifts need to be updated
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k_shift > nShifts
        k_shift = 1;
        tsh     = tic;

        [eqn, opts, oper, nShifts] = mess_lrradi_get_shifts(eqn, opts, oper, ...
                                                            W, Z, Y);
        
        timesh     = toc(tsh);
        out.timesh = out.timesh + timesh;

        if opts.shifts.info % print shifts
            fprintf('RADI Shifts:\n');
            disp(opts.shifts.p);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get current shift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pc       = opts.shifts.p( k_shift );
    out.p(k) = pc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [V, eqn, opts, oper] = mess_solve_shifted_system(eqn, opts, oper, pc, W);
    
    if opts.LDL_T
        V = V * diag(eqn.S_diag);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update low rank solution factor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isreal(pc)
        % The shift pc is real. Only perform a single step of the method.
        V     = real(V);
        VtB   = V' * eqn.BB;
        
        if opts.LDL_T
            Y_new = diag(eqn.S_diag);
        else
            Y_new = eye(p);
        end
        Y_new = Y_new + VtB * VtB';
        
        if opts.radi.compute_sol_fac || opts.radi.compute_sol_facpart
            % Only store part of Z used for shift generation.
            nZ = size(Z, 2);
            % Expand the Z matrix.
            if opts.radi.compute_sol_facpart
                ind = max(nZ-maxcolZ+p, 0)+1 : max(min(maxcolZ, nZ), 0);
                Z(:, 1:min(maxcolZ, k*p + nZ0)) = [Z(:, ind), sqrt(-2.0*pc)*V];
            else
                Z(:, nZ+1:nZ+p) = sqrt(-2.0*pc)*V;
            end
            
            if opts.radi.compute_sol_fac, Y = blkdiag(Y, Y_new); end
        end
        
        % Update the low-rank residual factor W.
        VY_newi = -2.0 * pc * (V / Y_new);
        if eqn.haveE
            VY_newi = oper.mul_E(eqn, opts, eqn.type, VY_newi, 'N');
        end
        W = W + VY_newi;
        
        % Update the K matrix.
        if eqn.type == 'T'            
            eqn.V(:, end-m+1:end) = eqn.V(:, end-m+1:end) + VY_newi * VtB;
        else
            eqn.U( : , end-m+1:end) = eqn.U( : , end-m+1:end) + VY_newi * VtB;
        end
    else
        % The shift pc is complex.
        % Perform a double step with the known solution for the conjugate
        % shift.
        V1 = sqrt(-2.0*real(pc)) * real(V);
        V2 = sqrt(-2.0*real(pc)) * imag(V);
        
        % Some auxiliary matrices.
        Vr = V1' * eqn.BB; 
        Vi = V2' * eqn.BB; 

        % Compute the new parts of low-rank approximate solution.
        sr = real(pc);
        si = imag(pc);
        sa = abs(pc);
        
        AA = [-sr/sa * Vr - si/sa * Vi; si/sa * Vr - sr/sa * Vi];
        BB = [Vr; Vi];
        CC = [si/sa * eye(p); sr/sa * eye(p)];
        
        if opts.LDL_T
            Y_new = blkdiag(diag(eqn.S_diag), 1/2 * diag(eqn.S_diag)) ...
                                             - 1/(4 * sr) * (AA * AA') ...
                                             - 1/(4 * sr) * (BB * BB') ...
                                             - 1/2 * (CC * (diag(eqn.S_diag) * CC'));
        else
            Y_new = blkdiag(eye(p), 1/2 * eye(p)) ...
                                   - 1/(4 * sr) * (AA * AA') ...
                                   - 1/(4 * sr) * (BB * BB') ...
                                   - 1/2 * (CC * CC');
        end
        
        if opts.radi.compute_sol_fac || opts.radi.compute_sol_facpart
            % Only store part of Z used for shift generation.
            nZ = size(Z, 2);
            % Expand the Z matrix. Z = [Z, V1, V2];
            if opts.radi.compute_sol_facpart
                ind = max(nZ-maxcolZ+2*p, 0) + 1:max(min(maxcolZ, nZ), 0);
                Z(:, 1:min( maxcolZ, (k+1)*p + nZ0)) = [Z(:, ind), V1, V2];
            else
                Z(:, nZ+1:nZ+2*p) = [V1, V2];
            end
            
            if opts.radi.compute_sol_fac
                Y = blkdiag(Y, Y_new);
            end
        end
        
        % Update the low-rank residual factor W.
        VY_newi = [V1, V2] / Y_new;
        if eqn.haveE
            VY_newi = oper.mul_E(eqn, opts, eqn.type, VY_newi, 'N');
        end
        
        W = W + sqrt(-2.0 * sr) * VY_newi(:, 1:p);
        
        % Update the K matrix.
        if eqn.type == 'T'
            eqn.V(:, end-m+1:end) = eqn.V(:, end-m+1:end) + VY_newi * [Vr; Vi];
        else
            eqn.U(:, end-m+1:end) = eqn.U(:, end-m+1:end) + VY_newi * [Vr; Vi];
        end
        
        % Forward the indices in the loop.
        k       = k + 1;
        k_shift = k_shift + 1;
        
        if k > 2
            res(k - 1) = res(k - 2);
        else
            res(k - 1) = 1;
        end
        out.p(k) = conj(pc);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.radi.res_tol
        % Low-rank residual norm computation.
        if opts.LDL_T
            res(k) = riccati_LR(W, [], opts, diag(eqn.S_diag), []) / res0;
        else
            res(k) = riccati_LR(W, [], opts, [], []) / res0;
        end
    end
    
    if opts.radi.rel_diff_tol
        if isreal(pc)
            nrmV = -2.0 * pc * sum(sum(V .^ 2));
        else
            % Complex double step means 2 blocks added.
            nrmV = sum(sum([V1, V2].^2 ));
        end
        nrmZ  = nrmZ + nrmV;
        rc(k) = sqrt(nrmV / nrmZ);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.radi.info
        if opts.radi.rel_diff_tol && opts.radi.res_tol
            fprintf(1, ...
                    ['RADI step: %4d pc: %e + %ei normalized residual: ' ...
                     '%e relative change in Z: %e\n'], ...
                    k, real(pc), imag(pc), res(k), rc(k));
        elseif opts.radi.res_tol
            fprintf(1, 'RADI step: %4d normalized residual: %e\n', ...
                    k, res(k));
        elseif opts.radi.rel_diff_tol
            fprintf(1, 'RADI step: %4d relative change in Z: %e\n', ...
                    k, rc(k));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if res(k) < opts.radi.res_tol
        break;
    end

    k       = k + 1;
    k_shift = k_shift + 1;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare output arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.niter = k - (k > opts.radi.maxiter);

if opts.radi.res_tol
    out.res = res(1:out.niter);
end

% warn the user if we have stopped before reaching the desired accuracy.
if (out.niter == opts.radi.maxiter) && ...
   not(out.res(end) < opts.radi.res_tol)

    warning('MESS:RADI:convergence',...
            ['LR-RADI was stopped by the maximum iteration count.',...
             ' Results may be inaccurate.'] );
end

if opts.radi.rel_diff_tol
    out.rc = rc(1:out.niter);
end

out.res_fact = W;

if eqn.type == 'T'
    out.K = eqn.V(:, end-m+1:end)';
else
    out.K = eqn.U(:, end-m+1:end)';
end

if opts.radi.compute_sol_fac && not(onlyK)

    if opts.radi.get_ZZt && not(opts.LDL_T)
        R     = chol(Y);
        out.Z = mess_column_compression(Z / R, 'N', [], opts.radi.trunc_tol, ...
                                        opts.radi.trunc_info);
    elseif opts.LDL_T
        Yinv           = Y \ eye(size(Y, 1));
        Yinv           = 0.5 * (Yinv + Yinv');
        [out.Z, out.D] = mess_column_compression(Z, 'N', Yinv, ...
                                                 opts.radi.trunc_tol, ...
                                                 opts.radi.trunc_info);
        out.Y          = out.D \ eye(size(out.D, 1));
    else
        out.Z = Z;
        out.Y = Y;
    end
end

out.res0 = res0;


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

% Delete short cuts for right hand-side and quadratic term.
eqn = rmfield(eqn, 'BB');
eqn = rmfield(eqn, 'CC');

if isfield(eqn, 'S_diag')
    eqn = rmfield(eqn, 'S_diag');
end

if isfield(opts.shifts, 'tmp')
    opts.shifts = rmfield(opts.shifts, 'tmp');
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalize required usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eqn.haveE
    [eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
end

[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.init_res_post(eqn, opts, oper);
