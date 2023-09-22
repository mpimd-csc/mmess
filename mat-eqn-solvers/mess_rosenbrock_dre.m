function [out, eqn, opts, oper] = mess_rosenbrock_dre(eqn, opts, oper)
%% [out, eqn, opts, oper] = mess_rosenbrock_dre(eqn, opts, oper)
%   LDL^T factored Rosenbrock method solving differential Riccati equation
%   E*d/dt X(t)*E' =-B*B' - E*X(t)*A' - A*X(t)*E' + E*X(t)*C'*C*X(t)*E' (N)
%   E'*d/dt X(t)*E =-C'*C - E'*X(t)*A - A'*X(t)*E + E'*X(t)*B*B'*X(t)*E (T)
%   backward in time.
%
%
% Input & Output
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation
%                       with A and E
%
% Output
%   out     struct contains solutions for every time step
%
% Input fields in struct eqn:
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   eqn.L0      dense (n x m4) matrix, initial L
%
%   eqn.D0      dense (m4 x m4) matrix, D in initial LDL^T factorization
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional)
%
%   eqn.haveE   possible  values: false, true
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
%   opts.rosenbrock.time_steps  possible  values: (N x 1) array
%                               array containing the time steps
%
%   opts.rosenbrock.stage       possible  values: 1, 2
%                               use 1-stage or 2-stage Rosenbrock method
%                               (optional, default: 1)
%
%   opts.rosenbrock.info        possible  values: 0, 1
%                               turn on (1) or off (0) the status output in
%                               every Rosenbrock iteration step
%                               (optional, default: 0)
%
%   opts.rosenbrock.gamma       possible  values: scalar
%                               Rosenbrock method coefficient
%                               (optional, default: 1+1/sqrt(2))
%
%   opts.rosenbrock.trunc_tol   possible  values: scalar > 0
%                               tolerance for LDL_T column compression
%                               (optional, default: eps*n)
%
%   opts.rosenbrock.trunc_info  possible values: 0, 1
%                               verbose mode for column compression
%                               (optional, default: 0)
%
%
%   opts.rosenbrock.
%   save_solution               possible  values: false, true
%                               save only K (0) or also the solution
%                               factors L and D (1)
%                               (optional, default: false)
%
%   opts.shifts.period          possible  values: integer > 0
%                               number of time steps that should pass
%                               until new shifts are computed
%                               (optional, default: 1)
%
% For the following fields see mess_lradi:
%
%   opts.adi.accumulateK
%
% Output fields in struct out:
%
%   out.Ks  cell array with matrix K for every time step
%
%   out.Ls  cell array with solution factor L for every time step
%           (only if opts.rosenbrock.save_solution = true)
%
%   out.Ds  cell array with solution factor D for every time step
%           (only if opts.rosenbrock.save_solution = true)
%
% If optional input arguments are missing they may be set to default values
% and a 'MESS:control_data' warning is printed. to turn warnings off use
% warning('OFF', 'MESS:control_data')
%
%
%   Note: Uses mess_lradi to solve inner Lyapunov Equations.
%
%   See also mess_lradi, mess_para, operatormanager.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for Rosenbrock Control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isfield(opts, 'rosenbrock')) || not(isstruct(opts.rosenbrock))
    mess_err(opts, 'control_data', ...
             'Rosenbrock control structure opts.rosenbrock missing.');
end % Single fields are checked below or inside mess_lradi

if not(isfield(opts.rosenbrock, 'time_steps'))
    mess_err(opts, 'control_data', 'opts.rosenbrock.time_steps is missing.');
end

if not(isfield(opts.rosenbrock, 'stage'))
    opts.rosenbrock.stage = 1;
end

if not(opts.rosenbrock.stage == 1) && not(opts.rosenbrock.stage == 2)
    mess_err(opts, 'control_data', ...
             ['opts.rosenbrock.stage has to be 1 or 2.', ...
              ' Other stages are not implemented']);
end

if not(isfield(opts.rosenbrock, 'info'))
    opts.rosenbrock.info = 0;
end

if not(isfield(opts.rosenbrock, 'trunc_tol'))
    opts.rosenbrock.trunc_tol = eps * oper.size(eqn, opts);
end

if not(isfield(opts.rosenbrock, 'trunc_info'))
    opts.rosenbrock.trunc_info = 0;
end

if not(isfield(opts.rosenbrock, 'gamma'))
    opts.rosenbrock.gamma = 1.0 + 1.0 / sqrt(2.0);
end

if not(isfield(opts.rosenbrock, 'save_solution'))
    opts.rosenbrock.save_solution = false;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for ADI control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts, 'adi')) || not(isstruct(opts.adi))
    mess_err(opts, 'control_data', 'ADI control structure opts.adi missing.');
end

if not(isfield(opts.adi, 'compute_sol_fac')) || ...
   not(islogical(opts.adi.compute_sol_fac)) || ...
   not(opts.adi.compute_sol_fac)

    mess_warn(opts, 'control_data', ...
              'Missing or Corrupted compute_sol_fac field. Switching to default.');
    opts.adi.compute_sol_fac = true;
end

if not(isfield(opts.adi, 'accumulateK')) || not(islogical(opts.adi.accumulateK))
    mess_warn(opts, 'control_data', ...
              'Missing or Corrupted accumulateK field. Switching to default: 1');
    opts.adi.accumulateK = true;
end

if not(isfield(opts.shifts, 'period'))
    opts.shifts.period = 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for shift computation control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts, 'shifts')) || not(isstruct(opts.shifts))
    mess_err(opts, 'control_data', 'shifts control structure opts.shifts missing.');
end

if not(isfield(opts.shifts, 'implicitVtAV'))
    opts.shifts.implicitVtAV = true;
end

if not(islogical(opts.shifts.implicitVtAV))
    mess_warn(opts, 'implicitVtAV', ...
              ['Missing or Corrupted implicitVtAV field.' ...
               'Switching to default (true).']);
    opts.shifts.implicitVtAV = true;
end

if not(opts.shifts.implicitVtAV)
    mess_warn(opts, 'implicitVtAV', ...
              ['implicitVtAV must be true for mess_rosenbrock_dre.', ...
               ' Switching to default (true).']);
    opts.shifts.implicitVtAV = true;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'LTV'))
    eqn.LTV = false;
end

if eqn.LTV
    mess_err(opts, 'control_data', ...
             ['non-autonomous differential Riccati', ...
              ' equation (eqn.LTV=1) is not supported']);
end

if not(isfield(eqn, 'C')) || not(isnumeric(eqn.C))
    mess_err(opts, 'control_data', ...
             'eqn.C is not defined or corrupted');
end

if not(isfield(eqn, 'B')) || not(isnumeric(eqn.B))
    mess_err(opts, 'control_data', 'eqn.B is not defined or corrupted');
end

if not(isfield(eqn, 'L0')) || not(isnumeric(eqn.L0))
    mess_warn(opts, 'control_data', ...
              ['Initial condition factor L0 is not defined or ', ...
               'corrupted. Setting it to the zero vector.']);
    n = oper.size(eqn);
    eqn.L0 = zeros(n, 1);
end

if not(isfield(eqn, 'D0')) || not(isnumeric(eqn.D0))
    mess_warn(opts, 'control_data', ...
              ['Initial condition factor D0 is not defined or ', ...
               'corrupted. Setting it to the identity matrix.']);
    eqn.D0 = eye(size(eqn.L0, 2));
end

if not(isfield(eqn, 'type'))
    eqn.type = 'N';
    mess_warn(opts, 'control_data', ...
              ['Unable to determine type of equation. ', ...
               'Falling back to type ''N''']);
elseif not(eqn.type == 'N') && not(eqn.type == 'T')
    mess_err(opts, 'equation_type', ...
             'Equation type must be either ''T'' or ''N''');
end

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');

if not(result)
    mess_err(opts, 'control_data', ...
             'system data is not completely defined or corrupted');
end

if eqn.type == 'T'
    q = size(eqn.C, 1);
else
    q = size(eqn.B, 2);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All checks done. Here comes the real work!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flip time_steps to run backwards in time
times = opts.rosenbrock.time_steps(end:-1:1);

L = eqn.L0;
D = eye(size(L, 2));

if eqn.type == 'T'
    if eqn.haveE
        K = oper.mul_E(eqn, opts, eqn.type, L * (L' * eqn.B), 'N');
    else
        K = L * (L' * eqn.B);
    end
else
    if eqn.haveE
        K = oper.mul_E(eqn, opts, eqn.type, L * (eqn.C * L)', 'N');
    else
        K = L * (eqn.C * L)';
    end
end

Iq = eye(q);

if opts.rosenbrock.save_solution
    out.Ls = {L}; % L of step 0
    out.Ds = {D}; % D of step 0
end

if eqn.type == 'T'
    eqn.V = K;
else
    eqn.U = K;
end

eqn.haveUV = true;
out.Ks = {K'}; % K of step 0

if eqn.type == 'T'
    eqn.W = eqn.C';
else
    eqn.W = eqn.B;
end

W = eqn.W;
opts.LDL_T = true;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 2:length(times)

    t = times(k);
    opts.rosenbrock.tau = times(k - 1) - times(k);

    if opts.rosenbrock.stage == 1
        if eqn.type == 'T'
            eqn.U = -eqn.B;
        else
            eqn.V = -eqn.C';
        end
    else
        if eqn.type == 'T'
            eqn.U = (-opts.rosenbrock.tau * opts.rosenbrock.gamma) * eqn.B;
        else
            eqn.V = (-opts.rosenbrock.tau * opts.rosenbrock.gamma) * eqn.C';
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update E^T * L
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if eqn.haveE
        EL = oper.mul_E(eqn, opts, eqn.type, L, 'N');
    else
        EL = L;
    end
    m = size(EL, 2);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update DLB block for S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if eqn.type == 'T'
        BLD_tmp = (eqn.B' * L) * D;
    else
        BLD_tmp = (eqn.C * L) * D;
    end

    if opts.rosenbrock.stage == 1
        BLD = BLD_tmp' * BLD_tmp + 1.0 / opts.rosenbrock.tau * D;

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update W and T
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        eqn.W = [W, EL];
    else % stage 2
        BLD = [zeros(m, m), D; ...
               D, -(BLD_tmp' * BLD_tmp)];

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update W and T
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        eqn.W =  [W, oper.mul_A(eqn, opts, eqn.type, L, 'N'), EL];
    end

    eqn.T = blkdiag(Iq, BLD);
    [eqn.W, eqn.T] = mess_column_compression(eqn.W, 'N', eqn.T, ...
                                             opts.rosenbrock.trunc_tol, ...
                                             opts.rosenbrock.trunc_info);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute new ADI shifts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if not(mod(k - 2, opts.shifts.period))
        opts.shifts.p = mess_para(eqn, opts, oper);
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.rosenbrock.stage == 1
        adiout = mess_lradi(eqn, opts, oper);

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform column compression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [L, D] = mess_column_compression(adiout.Z, 'N', adiout.D, ...
                                         opts.rosenbrock.trunc_tol, ...
                                         opts.rosenbrock.trunc_info);

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract or construct K
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if opts.adi.accumulateK
            K = adiout.Knew;
        else
            if eqn.type == 'T'
                if eqn.haveE
                    K = oper.mul_E(eqn, opts, eqn.type, L, 'N') * ...
                        (D * (L' * eqn.B));
                else
                    K = L * (D * (L' * eqn.B));
                end

            else
                if eqn.haveE
                    K = oper.mul_E(eqn, opts, eqn.type, L, 'N') * ...
                        (D * (L' * eqn.C'));
                else
                    K = L * (D * (L' * eqn.C'));
                end
            end
        end
    else % stage = 2
        adiout1 = mess_lradi(eqn, opts, oper);

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform column compression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [T1, D1] = mess_column_compression(adiout1.Z, 'N', adiout1.D, ...
                                           opts.rosenbrock.trunc_tol, ...
                                           opts.rosenbrock.trunc_info);

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update RHS for second equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if eqn.haveE
            ET1 = oper.mul_E(eqn, opts, eqn.type, T1, 'N');
        else
            ET1 = T1;
        end
        eqn.W = ET1;

        if eqn.type == 'T'
            BT1E_tmp = (eqn.B' * T1) * D1;
        else
            BT1E_tmp = (eqn.C * T1) * D1;
        end

        BT1D = opts.rosenbrock.tau * opts.rosenbrock.tau * ...
               (BT1E_tmp' * BT1E_tmp) + ...
               (2.0 - 1.0 / opts.rosenbrock.gamma) * D1;

        eqn.T = BT1D;

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute new ADI shifts
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if not(mod(k - 2, opts.shifts.period))
            opts.shifts.p = mess_para(eqn, opts, oper);
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % solve second equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        adiout2 = mess_lradi(eqn, opts, oper);

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct new X = L * D * L^T
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L = [L, T1, adiout2.Z]; %#ok<AGROW>
        D = blkdiag(D, opts.rosenbrock.tau * ...
                    (2.0 - 1.0 / (2.0 * opts.rosenbrock.gamma)) * D1, ...
                    -opts.rosenbrock.tau / 2.0 * adiout2.D);

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform column compression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [L, D] = mess_column_compression(L, 'N', D, ...
                                         opts.rosenbrock.trunc_tol, ...
                                         opts.rosenbrock.trunc_info);
        if eqn.type == 'T'
            eqn.W = eqn.C';
        else
            eqn.W = eqn.B;
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct new K
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if opts.adi.accumulateK
            K = K + opts.rosenbrock.tau * ...
                (2.0 - 1.0 / (2.0 * opts.rosenbrock.gamma)) * ...
                adiout1.Knew - opts.rosenbrock.tau / 2.0 * adiout2.Knew;
        else
            if eqn.type == 'T'
                if eqn.haveE
                    K = oper.mul_E(eqn, opts, eqn.type, L, 'N') * ...
                        (D * (L' * eqn.B));
                else
                    K = L * (D * (L' * eqn.B));
                end
            else
                if eqn.haveE
                    K = oper.mul_E(eqn, opts, eqn.type, L, 'N') * ...
                        (D * (L' * eqn.C'));
                else
                    K = L * (D * (L' * eqn.C'));
                end
            end
        end
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opts.rosenbrock.info
        mess_fprintf(opts, 'Rosenbrock step %4d s\n', t);
    end

    if opts.rosenbrock.save_solution
        out.Ls = [{L}, out.Ls];
        out.Ds = [{D}, out.Ds];
    end

    out.Ks = [{K'}, out.Ks];

    if eqn.type == 'T'
        eqn.V = K;
    else
        eqn.U = K;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
