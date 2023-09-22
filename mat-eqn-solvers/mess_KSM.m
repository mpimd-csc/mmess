function [out, eqn, opts, oper] = mess_KSM(eqn, opts, oper)
% function [out, eqn, opts, oper] = mess_KSM(eqn, opts, oper)
% Solve continuous-time Lyapunov equations with sparse coefficients
%
%   eqn.type = 'N' -> A*Z*Z'*E' + E*Z*Z'*A' + B*B' = 0 (N)
%   eqn.type = 'T' -> A'*Z*Z'*E + E'*Z*Z'*A + C'*C = 0 (T)
%
% or continuous time algebraic Riccati equation
%
%   eqn.type = 'N' -> A*X*E' + E*X*A' - E*X*C'*C*X*E' + B*B' = 0 (N)
%   eqn.type = 'T' -> A'*X*E + E'*X*A - E'*X*B*B'*X*E + C'*C = 0 (T)
%
% by th extended and rational Krylov subspace method.
%
% Input & Output
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operations
%                       with A and E
%
% Output
%   out     struct contains the low-rank factor of the solution (and other
%           info)
%
% Input fields in struct eqn:
%
%   eqn.B       dense (n x q) matrix B
%
%   eqn.C       dense (p x n) matrix C
%
%
%   eqn.type    possible values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional, default 'N')
%
%   eqn.haveE   possible values: false, true
%               if haveE = false: matrix E is assumed to be the identity
%               (optional)
%
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices. For the second order types these are given
%   implicitly by the M, D, K matrices stored in eqn.M_, eqn.E_ and eqn.K_,
%   respectively.
%
% Input fields in struct opts:
%   opts.KSM.space       possible values:
%                          'EK': the extended Krylov subspace method is applied
%                          'RK': the rational Krylov subspace method is applied
%
%   opts.KSM.type_eq     possible values:
%                          "LE": if we are solving a Lyapunov
%                                equation
%                          "CARE": if we are solving an algebraic
%                                 Riccati equation
%
%   opts.KSM.symmetric   possible values
%                          1: if A is symmetric
%                             any other value if A is nonsymmetric
%                             REMARK: full orthogonalization is always
%                             performed to improve stability
%
%   opts.KSM.maxiter      maximum number of iteration allowed
%
%   opts.KSM.res_tol      threshold for the relative residual norm.
%                         If res < opts.tol * norm(rhs) we stop the algorithm
%
%   opts.KSM.trunc_tol     possible values: scalar > 0
%                          tolerance for rank truncation of the
%                          low-rank solutions (aka column compression)
%                          (optional, default: eps*n)
%
%   opts.KSM.trunc_info    possible values: 0, 1
%                          verbose mode for column compression
%                          (optional, default: 0)
%
% Output fields in struct out:
%
%   out.Z, out.D     the approximate solution is represented as
%                                out.Z*out.D*out.Z'
%                    where out.Z is n-by-t, t<<n,
%                    and out.D is a diagonal t-by-t matrix
%   out.niter        number of iterations
%   out.res          array of relative RADI residual norms
%
% If optional input arguments are missing they may be set to default values
% and a 'MESS:control_data' warning is printed. to turn warnings off use
% warning('OFF', 'MESS:control_data')
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check field opts.KSM
if not(isfield(opts, 'KSM')) || not(isstruct(opts.KSM))
    mess_err(opts, 'control_data', ...
             'No KSM control data found in options structure.');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check info parameter for output verbosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.KSM, 'info'))
    opts.KSM.info = 0;
else
    if not(isnumeric(opts.KSM.info)) && not(islogical(opts.KSM.info))
        mess_err(opts, 'info', ...
                 'opts.KSM.info parameter must be logical or numeric.');
    end
end

if not(isfield(opts, 'shifts')) || not(isfield(opts.shifts, 'info'))
    opts.shifts.info = 0;
else
    if not(isnumeric(opts.shifts.info)) && not(islogical(opts.shifts.info))
        mess_err(opts, 'info', ...
                 'opts.shifts.info parameter must be logical or numeric.');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check stopping parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.KSM, 'maxiter')) || not(isnumeric(opts.KSM.maxiter))
    mess_warn(opts, 'control_data', ...
              ['Missing or Corrupted opts.KSM.maxiter field. ', ...
               'Switching to default: 100']);
    opts.KSM.maxiter = 100;
end

if not(isfield(opts.KSM, 'res_tol')) || not(isnumeric(opts.KSM.res_tol))
    mess_warn(opts, 'control_data', ...
              ['Missing or Corrupted opts.KSM.res_tol field. ', ...
               'Switching to default: 1e-6']);
    opts.KSM.res_tol = 1e-6;
end

if not(isfield(opts.KSM, 'trunc_tol'))
    opts.KMS.trunc_tol = eps * oper.size(eqn, opts);
end
if not(isfield(opts.KSM, 'trunc_info'))
    opts.KSM.trunc_info = 0;
end

if not(isfield(opts.KSM, 'comp_res')) || not(isnumeric(opts.KSM.comp_res))
    mess_warn(opts, 'control_data', ...
              ['Missing or Corrupted opts.KSM.comp_res field. ', ...
               'Switching to default: 10']);
    opts.KSM.comp_res = 10;
end

if not(isfield(opts, 'norm')) || ...
        (not(strcmp(opts.norm, 'fro')) && ...
         (not(isnumeric(opts.norm)) || not(opts.norm == 2)))
    mess_warn(opts, 'control_data', ...
              ['Missing or Corrupted opts.norm field. ', ...
               'Switching to default: ''fro''']);
    opts.norm = 'fro';
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'type'))
    eqn.type = 'N';
    mess_warn(opts, 'control_data', ...
              ['Unable to determine type of equation.'...
               'Falling back to type ''N''']);
elseif not(eqn.type == 'N') && not(eqn.type == 'T')
    mess_err(opts, 'equation_type', ...
             'Equation type must be either ''T'' or ''N''');
end

% set flag 0 if E does not exist
if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
    mess_warn(opts, 'control_data', ...
              ['Missing or Corrupted eqn.haveE field.', ...
               'Switching to default: 0']);
end

[erg, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');

if not(erg)
    mess_err(opts, 'control_data', ...
             'system data is not completely defined or corrupted');
end

% initialize solving and multiplying with A for EKSM and shift computation
% in RKSM, shifted solves for RKSM are handled inside the space expansion
% routine, as this needs to be done per shift.
[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_A_pre(eqn, opts, oper);

% statespace transformation is required in any case:
[eqn, opts, oper] = oper.dss_to_ss_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.ss_to_dss_pre(eqn, opts, oper);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check initial shifts and the computation of complex/real shifts for RKSM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(opts.KSM.space, 'RK')
    if not(isfield(opts.KSM, 'init_shifts'))
        opts_eigs.tol = 1e-1;
        A = @(b)(oper.mul_A(eqn, opts, 'N', b, 'N'));
        n = oper.size(eqn, opts);
        try % new eigs syntax if available (R2017b and newer)
            opts.KSM.init_shifts(1) = ...
                -eigs(A, n, 1, 'largestabs', opts_eigs);
            opts.KSM.init_shifts(2) = ...
                -eigs(A, n, 1, 'smallestabs', opts_eigs);
        catch % old syntax otherwise (octave and pre-R2017b)
            opts.KSM.init_shifts(1) = -eigs(A, n, 1, 'lm', opts_eigs);
            opts.KSM.init_shifts(2) = -eigs(A, n, 1, 'sm', opts_eigs);
        end
    end

    if not(isfield(opts.KSM, 'type_shifts'))
        if isfield(opts.KSM, 'symmetric') && opts.KSM.symmetric
            opts.KSM.type_shifts = 'real';
            mess_warn(opts, 'control_data', ...
                      ['Missing or Corrupted opts.KSM.type_shifts ', ...
                       'field. Switching to default: ''real''']);
        else
            opts.KSM.type_shifts = 'complex';
            mess_warn(opts, 'control_data', ...
                      ['Missing or Corrupted opts.KSM.type_shifts ', ...
                       'field. Switching to default: ''complex''']);
        end
    end

    % If CARE, we compute the shifts by using the eigenvalues of the
    % projected closed-loop operator and not the Ritz values
    if strcmp(opts.KSM.type_eqn, 'CARE')
        if not(isfield(opts.KSM, 'CARE_shifts'))
            opts.KSM.CARE_shifts = 'Ritz_closedloop';
            mess_warn(opts, 'control_data', ...
                      ['Missing or Corrupted opts.KSM.CARE_shifts ', ...
                       'field. Switching to default: ''Ritz_closedloop''']);
        end
    end

end

%%
% check whether LDL^T formulation should be used
if not(isfield(opts, 'LDL_T'))
    opts.LDL_T = false;
end
% check for or set proper right hand side in eqn.W.
% We use eqn.W as starting block
if opts.LDL_T
    % RHS has form W* T * W'
    if not(isfield(eqn, 'W')) || not(isnumeric(eqn.W))
        mess_err(opts, 'control_data', ...
                 'eqn.W is not defined or corrupted');
    end
    if not(isfield(eqn, 'T')) || not(isnumeric(eqn.T))
        mess_err(opts, 'control_data', ...
                 'eqn.T is not defined or corrupted');
    end

else
    if not(isfield(eqn, 'W'))
        if eqn.type == 'N'
            eqn.W = eqn.B;
        else
            eqn.W = eqn.C';
        end
        eqn.T = eye(size(eqn.W, 2));
    end
end

% Compute B'*B (small matrix) as \|B*B'\|_F=\|B'*B\|_F for residual
% normalization. Note that the residual is computed for the original
% non transformed equation. So we use the original B and C
norm_rhs = sqrt(trace(eqn.T' * (eqn.W' * eqn.W) * ...
                      eqn.T * (eqn.W' * eqn.W)));
p = size(eqn.W, 2);
tol = opts.KSM.res_tol * norm_rhs;

% state space transformation (if needed)
if eqn.haveE
    eqn.ssW = full(oper.dss_to_ss(eqn, opts, 'L', 'N', eqn.W, 'N'));
    if strcmp(opts.KSM.type_eqn, 'CARE')
        if eqn.type == 'N'
            eqn.ssC = full(oper.dss_to_ss(eqn, opts, 'U', 'T', ...
                                          eqn.C, 'T'))';
        else
            eqn.ssB = full(oper.dss_to_ss(eqn, opts, 'L', 'N', ...
                                          eqn.B, 'N'));
        end
    end
else
    eqn.ssW = eqn.W;
    eqn.ssB = eqn.B;
    eqn.ssC = eqn.C;
end

% Orthonormalize the first basis block [B,A^{-1}*B]
opts.KSM.compute_struct.V = [];
[eqn, opts, oper] = KSM_space_expansion(eqn, opts, oper);

for j = 1:opts.KSM.maxiter

    % Expand the space
    [eqn, opts, oper] = KSM_space_expansion(eqn, opts, oper);

    if not(mod(j, opts.KSM.comp_res))
        % Solve the projected problem
        [eqn, opts, oper] = mess_solve_projected_eqn(eqn, opts, oper, ...
                                                     'KSM', ...
                                                     opts.KSM.type_eqn);

        % Compute the residual norm
        [eqn, opts, oper] = KSM_compute_res(eqn, opts, oper);

        if opts.KSM.info
            mess_fprintf(opts, ['KSM residual at iteration %3d is:  ', ...
                                '%10e (absolute)  %10e (normalized)\n'], ...
                         j, opts.KSM.compute_struct.res(end), ...
                         opts.KSM.compute_struct.res(end) / norm_rhs);
        end

        if opts.KSM.compute_struct.res(end) < tol
            break
        end
    end

end

% warn the user if we have stopped before reaching the desired accuracy.
if j == opts.KSM.maxiter
    mess_warn(opts, 'convergence', ...
              ['KSM reached maximum iteration number. ', ...
               'Results may be inaccurate!']);
end

% compute the eigen decomposition of Y
[L, D] = eig(opts.KSM.compute_struct.Y);
% and perform a low-rank truncation
if opts.KSM.trunc_info
    [L, D] = mess_column_compression(L, 'N', D, opts.KSM.trunc_tol);
end

% The solution is represented as Z*D*Z', D diagonal
if strcmp('EK', opts.KSM.space)
    out.Z = opts.KSM.compute_struct.V(:, 1:end - 2 * p) * L;
elseif strcmp('RK', opts.KSM.space)
    out.Z = opts.KSM.compute_struct.V(:, 1:end - p) * L;
end

if eqn.haveE
    switch eqn.type
        case 'N'
            out.Z = oper.dss_to_ss(eqn, opts, 'U', 'N', out.Z, 'N');
        case 'T'
            out.Z = oper.dss_to_ss(eqn, opts, 'L', 'T', out.Z, 'N');
    end
end
out.D = D;
out.niter = j;
out.res = opts.KSM.compute_struct.res / norm_rhs;

[eqn, opts, oper] = oper.dss_to_ss_post(eqn, opts, oper);
[eqn, opts, oper] = oper.ss_to_dss_post(eqn, opts, oper);

if strcmp('EK', opts.KSM.space)
    [eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
    [eqn, opts, oper] = oper.sol_A_pre(eqn, opts, oper);
end
