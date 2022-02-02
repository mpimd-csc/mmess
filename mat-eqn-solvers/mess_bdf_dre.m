function [out, eqn, opts, oper] = mess_bdf_dre(eqn, opts, oper)
%% function [out, eqn, opts, oper] = mess_bdf_dre(eqn, opts, oper)
%   LDL^T factored BDF method solving differential Riccati equation
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
%   eqn.B       dense (n x m1) matrix B (if eqn.LTV = 0)
%               (optional)
%
%   eqn.B_time  function handle with scalar input t returning a dense
%               (n x m1) matrix B (if eqn.LTV = 1)
%               (optional)
%
%   eqn.C       dense (m2 x n) matrix C (if eqn.LTV = 0)
%               (optional)
%
%   eqn.C_time  function handle with scalar input t returning a dense
%               (m2 x n) matrix C (if eqn.LTV = 1)
%               (optional)
%
%   eqn.L0      dense (n x m4) matrix, initial L
%               (optional)
%
%   eqn.D0      dense (m4 x m4) matrix, D in initial LDL^T factorization
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional)
%
%   eqn.haveE   possible  values: 0, 1, false, true
%               if haveE = 0: matrix E is assumed to be the identity
%               (optional)
%
%   eqn.LTV     possible  values: 0, 1, false, true
%               indicates autonomous (false) or
%               non-autonomous (true) differential
%               Riccati equation
%               (optional, default: 0)
%
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices. For the second order types these are given
%   implicitly by the M, D, K matrices stored in eqn.M_, eqn.E_ and eqn.K_,
%   respectively.
%
% Input fields in struct opts:
%   opts.bdf.time_steps             possible  values: (N x 1) array
%                                   array containing the time steps
%
%   opts.bdf.step                   possible  values: 1, 2, 3, 4
%                                   use p-step  BDF method with 1<=p<=4
%                                   (optional, default: 1)
%
%   opts.bdf.info                   possible  values: 0, 1, false, true
%                                   turn on (1) or off (0) the status
%                                   output in every BDF iteration step
%                                   (optional, default: 0)
%
%   opts.bdf.save_solution          possible  values: 0, 1, false, true
%                                   save only K (0) or also the solution
%                                   factors L and D (1)
%                                   (optional, default: 0)
%
%   opts.bdf.trunc_tol              possible values: scalar > 0
%                                   tolerance for LDL_T column compression
%                                   (optional, default: eps*n)
%
%   opts.bdf.trunc_info             possible values: 0, 1, false, true
%                                   verbose mode for column compression
%                                   (optional, default: 0)
%
%   opts.nm                         struct for control parameters of
%                                   mess_lrnm
%
%   opts.adi                        struct for control parameters of
%                                   mess_lradi
%
% Output fields in struct out:
%
%   out.Ks  cell array with matrix K for every time step
%
%   out.Ls  cell array with solution factor L for every time step
%           (only if opts.bdf.save_solution = 1)
%
%   out.Ds  cell array with solution factor D for every time step
%           (only if opts.bdf.save_solution = 1)
%
%
% If optional input arguments are missing they may be set to default values
% and a 'MESS:control_data' warning is printed. to turn warnings off use
% warning('OFF', 'MESS:control_data')
%
%   Note: uses mess_lrnm to solve inner Riccati Equations
%
%   See also mess_lrnm, mess_lradi, mess_para, operatormanager.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for BDF Control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isfield(opts,'bdf')) || not(isstruct(opts.bdf))
    error('MESS:control_data','BDF control structure opts.bdf missing.');
end % Single fields are checked below or inside mess_lradi

if not(isfield(opts.bdf,'time_steps')) || not(isnumeric(opts.bdf.time_steps)) ...
   || isscalar(opts.bdf.time_steps)
    error('MESS:control_data','opts.bdf.time_steps is missing.');
end

opts.t0 = opts.bdf.time_steps(1);
t1 = opts.bdf.time_steps(2);
tend = opts.bdf.time_steps(end);
opts.bdf.tau = t1 - opts.t0;
tau_original = opts.bdf.tau;
eq_err = norm(opts.bdf.time_steps - linspace(opts.t0, tend, length(opts.bdf.time_steps)));

if eq_err > (eps * length(opts.bdf.time_steps))
    error('MESS:control_data', ...
          'opts.bdf.time_steps has to contain equidistant time steps.');
end

if not(isfield(opts.bdf,'step')),  opts.bdf.step = 1; end

if rem(opts.bdf.step, 1) || (opts.bdf.step < 1) || (opts.bdf.step > 4)
    error('MESS:control_data','opts.bdf.step has an invalid value.');
end

if not(isfield(opts.bdf,'info')), opts.bdf.info = 0; end

if not(isfield(opts.bdf,'trunc_tol'))
    opts.bdf.trunc_tol = eps * oper.size(eqn, opts);
end

if not(isfield(opts.bdf, 'trunc_info')), opts.bdf.trunc_info = 0; end

if not(isfield(opts.bdf,'save_solution')), opts.bdf.save_solution = 0; end

if opts.bdf.step > 2
    if not(isfield(opts.bdf,'startup_iter') )
        opts.bdf.startup_iter = 8;
    end

    if rem(opts.bdf.startup_iter, 1)
        error('MESS:control_data',...
              'opts.bdf.startup_iter has an invalid value.');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for Newton control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts,'nm')) || not(isstruct(opts.nm))
    error('MESS:control_data','Newton control structure opts.nm missing.');
end

if isfield(opts.nm,'res_tol') && isnumeric(opts.nm.res_tol)
    out.res = [];
end

if isfield(opts.nm,'rel_diff_tol') && isnumeric(opts.nm.rel_diff_tol)
    out.rc = [];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for ADI control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts,'adi')) || not(isstruct(opts.adi))
    error('MESS:control_data','ADI control structure opts.adi missing.');
end

if not(isfield(opts.adi,'compute_sol_fac')) || ...
   not(isnumeric(opts.adi.compute_sol_fac)) || ...
   not(opts.adi.compute_sol_fac)

    warning('MESS:compute_sol_fac', ...
           'Missing or Corrupted compute_sol_fac field. Switching to default.');
    opts.adi.compute_sol_fac = 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for shift computation control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts,'shifts')) || not(isstruct(opts.shifts))
    error('MESS:control_data','shifts control structure opts.shifts missing.');
end

if not(isfield(opts.shifts,'implicitVtAV'))
    opts.shifts.implicitVtAV = true;
end

if not(isnumeric(opts.shifts.implicitVtAV)) && ...
   not(islogical(opts.shifts.implicitVtAV))
    warning('MESS:implicitVtAV', ...
       'Missing or Corrupted implicitVtAV field. Switching to default (true).');
    opts.shifts.implicitVtAV = true;
end

if not(opts.shifts.implicitVtAV)
    warning('MESS:implicitVtAV', ...
    'implicitVtAV must be true for mess_bdf_dre. Switching to default (true).');
    opts.shifts.implicitVtAV = true;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'LTV')), eqn.LTV = 0; end

if isfield(opts,'LDL_T') && (opts.LDL_T == 0)
    warning('MESS:control_data',['The BDF code does only support ' ...
                                 'LDL_T solutions.\n Setting opts.LDL_T=1']);
end

opts.LDL_T = 1;

if not(isfield(eqn,'type'))
    eqn.type='N';
    warning('MESS:control_data',['Unable to determine type of equation.'...
                                 'Falling back to type ''N''']);

elseif not(eqn.type=='N') && not(eqn.type=='T')
    error('MESS:equation_type','Equation type must be either ''T'' or ''N''');
end

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');

if eqn.LTV
    if isfield(oper, 'eval_matrix_functions')
        [eqn, opts, oper] = oper.eval_matrix_functions(eqn, opts, oper, ...
                                                       opts.bdf.time_steps(end));
    else
        error('MESS:missing_feature', ['The function eval_matrix_functions is ', ...
                                       'required for LTV problems, but it has ', ...
                                       'not yet been implemented for this set ', ...
                                       'of USFS functions']);
    end
end

if not(isfield(eqn, 'C')) || not(isnumeric(eqn.C))
    error('MESS:control_data', 'eqn.C is not defined or corrupted');
end

if not(isfield(eqn, 'B')) || not(isnumeric(eqn.B))
    error('MESS:control_data', 'eqn.B is not defined or corrupted');
end

if not(isfield(eqn, 'L0')) || not(isnumeric(eqn.L0))
    warning('MESS:control_data', ...
            ['Initial condition factor L0 is not defined or corrupted.',...
             ' Setting it to the zero vector.']);
    n = oper.size(eqn, opts);
    eqn.L0 = zeros(n,1);
end

if not(isfield(eqn, 'D0')) || not(isnumeric(eqn.D0))
    warning('MESS:control_data', ...
            ['Initial condition factor D0 is not defined or corrupted.',...
             ' Setting it to the identity matrix.']);
    eqn.D0 = eye(size(eqn.L0, 2));
end

if not(result)
    error('MESS:control_data','system data is not completely defined or corrupted');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize BDF parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = zeros(6, 6);
alpha(1, 1) = -1;
alpha(1 : 2, 2) = [-4. / 3.; 1. / 3.];
alpha(1 : 3, 3) = [-18. / 11.; 9. / 11.; -2. / 11.];
alpha(1 : 4, 4) = [-48. / 25.; 36. / 25.; -16. / 25.; 3. / 25.];
beta = [1; 2 / 3; 6 / 11; 12 / 25; 60 / 137; 60 / 147];
opts.bdf.beta = beta(opts.bdf.step);
step = opts.bdf.step;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

times = opts.bdf.time_steps(end:-1:1);
extra_steps = 0;

if step == 3
    % add small time steps to initialize L with correct order
    tau_small = opts.bdf.tau / 2^(opts.bdf.startup_iter);
    tau_old = tau_small;
    opts.bdf.tau = tau_old;
    times_extra = zeros(1, opts.bdf.startup_iter + 3);
    times_extra(1) = times(1); % copy endtime
    times_extra(2) = times(1) - tau_small; % one BDF 1 step
    times_extra(3:end) = times(1) - ...
                         tau_small * 2.^(1:opts.bdf.startup_iter + 1); % BDF 2 steps
    times = [times_extra, times(4 : end)];
    extra_steps = length(times_extra);

elseif step == 4
    % add small time steps to initialize L with correct order
    tau_small = opts.bdf.tau / 2^(opts.bdf.startup_iter + 1);
    tau_old = tau_small;
    opts.bdf.tau = tau_old;
    times_extra = zeros(1, 2 * opts.bdf.startup_iter + 5);
    times_extra(1) = times(1); % copy endtime
    times_extra(2) = times(1) - tau_small; % one BDF 1 step
    times_extra(3) = times(1) - 2 * tau_small; % one BDF 2 step
    % BDF 3 steps
    times_extra(4 : 2 : end) = times(1) - tau_small * ...
        (2.^(1:opts.bdf.startup_iter + 1) + 2.^(0:opts.bdf.startup_iter));
    times_extra(5 : 2 : end) = times(1) - tau_small * ...
        (2.^(2:opts.bdf.startup_iter + 2));
    times = [times_extra, times(1) - 3 * tau_original, times(5 : end)];
    extra_steps = length(times_extra) + 1;
end
L = eqn.L0;
D = eye(size(L, 2));

if eqn.type == 'T'
    K0 = oper.mul_E(eqn, opts, eqn.type, (L * (L' * eqn.B)), 'N');
else
    K0 = oper.mul_E(eqn, opts, eqn.type, (L * (L' * eqn.C')), 'N');
end

L_lengths = zeros(opts.bdf.step, 1);
Ds = {};
Iq = eye(q);

if opts.bdf.save_solution
    out.Ls = {L}; % L of step 0
    out.Ds = {D}; % D of step 0
elseif eqn.LTV % for LTV 2 step BDF save last L anyway
    out.Ls = {L}; % L of step 0
end

out.Ks = {K0'}; % K of step 0

% FIXME: dead code
%     end
%     out.Ks = zeros(p,oper.size(eqn)); % K of step 0
if eqn.type == 'T'
    C = eqn.C;
else
    B = eqn.B;
end

ETL = [];

if not(eqn.LTV)
    ETL_0 = oper.mul_E(eqn, opts, eqn.type, L, 'N');
end

Ds_0 = D;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 2 : length(times)
    t = times(k);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set order and beta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    opts.bdf.step = min(step, k - 1);
    tau_old = opts.bdf.tau;
    opts.bdf.tau = times(k - 1) - times(k);

    if (k <= extra_steps) && (opts.bdf.step == step)
        % additional smaller time steps with one order less
        opts.bdf.step = opts.bdf.step - 1;
    end

    opts.bdf.beta = beta(opts.bdf.step);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if eqn.LTV
        [eqn, opts, oper] = oper.eval_matrix_functions(eqn, opts, oper, t);
        if eqn.type == 'T'
            C = eqn.C;
        else
            B = eqn.B;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set order and beta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     opts.bdf.step = min(step, j - 1); % FIXME: dead code

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update E' * L
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if abs(tau_old - opts.bdf.tau) > (1e+4 * eps) % FIXME: add comment for magic number
        % in startup phase to initialize L with correct order, tau was
        % doubled from last time step to this, remove intermediate step
        if step == 3
            if eqn.LTV || opts.bdf.save_solution
                out.Ls = [out.Ls(1 : end - 2), out.Ls(end)];
                if opts.bdf.save_solution
                    out.Ds = [out.Ds(1 : end - 2), out.Ds(end)];
                end
            end

            if not(eqn.LTV)
                ETL( : , 1 : L_lengths(1)) = [];
            end

            L_lengths(1) = L_lengths(2);
            L_lengths(2) = 0;
            out.Ks = [out.Ks(1 : end - 2), out.Ks(end)];
            Ds = [Ds(1 : end - 2), Ds(end)];

            if isfield(out, 'res')
                out.res = [out.res(1 : end - 2), out.res(end)];
            end

            if isfield(out, 'rc')
                out.rc = [out.rc(1 : end - 2), out.rc(end)];
            end

        elseif step == 4
            if eqn.LTV || opts.bdf.save_solution
                out.Ls = [out.Ls(1 : end - 4), out.Ls(end - 2), out.Ls(end)];

                if opts.bdf.save_solution
                    out.Ds = [out.Ds(1 : end - 4), out.Ds(end - 2), out.Ds(end)];
                end
            end

            if not(eqn.LTV)
                ETL( : , sum(L_lengths(1:2)) + 1 : sum(L_lengths(1:3))) = [];
                ETL( : , 1 : L_lengths(1)) = [];
                ETL = [ETL, ETL_0]; %#ok<AGROW>
%                 L_lengths(4) = size(ETL_0, 2); % FIXME: dead code
            end

            Ds = [Ds, {Ds_0}]; %#ok<AGROW>

            L_lengths = [L_lengths([2, 4]); 0; 0];
            out.Ks = [out.Ks(1 : end - 4), out.Ks(end - 2), out.Ks(end)];
            Ds = [Ds(1 : end - 4), Ds(end - 2), Ds(end)];

            if isfield(out, 'res')
                out.res = [out.res(1 : end - 4), out.res(end - 2), out.res(end)];
            end

            if isfield(out, 'rc')
                out.rc = [out.rc(1 : end - 4), out.rc(end - 2), out.rc(end)];
            end
        end
% FIXME: dead code
%     elseif step == 4 && opts.bdf.step == step - 1 && j > 4 ...X
%             && not(eqn.LTV)
%         ETL_0 = ETL( : , end - L_lengths(opts.bdf.step) + 1 : end);
%         Ds_0 = Ds{end};
    end
    %     if not(isempty(L)) % FIXME: dead code
    if eqn.LTV

        for s = 2 : opts.bdf.step
            L = [L, out.Ls{s}]; %#ok<AGROW>
        end

        ETL = [];
    end

    if isempty(ETL)
        ETL = oper.mul_E(eqn, opts, eqn.type, L, 'N');
    else
        ETL = [oper.mul_E(eqn, opts, eqn.type, L, 'N'), ...
            ETL( : , 1 : end - L_lengths(opts.bdf.step))];
        L_lengths(2 : step) = L_lengths(1 : end - 1);
    end
    L_lengths(1) = size(L, 2);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update D blocks for S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    alphaDs = -alpha(1, opts.bdf.step) * D;

    for l = 2 : opts.bdf.step

        if size(Ds, 2) >= l - 1
            alphaDs = blkdiag(alphaDs, -alpha(l, opts.bdf.step) * Ds{l - 1});
        end
    end

    if size(Ds, 2) == opts.bdf.step
        Ds = {D, Ds{1 : end - 1}};
    else
        Ds = [{D}, Ds]; %#ok<AGROW>
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update C/B and S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if eqn.type == 'T'
        eqn.C = [C; ETL'];
    else
        eqn.B = [B, ETL];
    end

    eqn.S = blkdiag(opts.bdf.tau * opts.bdf.beta * Iq, alphaDs);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform column compression
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if eqn.type == 'T'
        [eqn.C, eqn.S] = mess_column_compression(eqn.C, 'T', eqn.S, ...
            opts.bdf.trunc_tol, opts.bdf.trunc_info);
    else
        [eqn.B, eqn.S] = mess_column_compression(eqn.B, 'N', eqn.S, ...
            opts.bdf.trunc_tol, opts.bdf.trunc_info);
    end

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform the actual step computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nmout = mess_lrnm(eqn, opts, oper);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform column compression
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [L, D] = mess_column_compression(nmout.Z, 'N', nmout.D, ...
                                     opts.bdf.trunc_tol, opts.bdf.trunc_info);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print status information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if opts.bdf.info
        fprintf('BDF step %4d s\n', t);
        fprintf('\t Newton steps: %2d \n', nmout.niter);
        fprintf('\t Rank %d\n', size(Ds{1},1));
    end

    if opts.bdf.save_solution
        out.Ls = [{L}, out.Ls];
        out.Ds = [{D}, out.Ds];
    elseif eqn.LTV % for LTV 2 step BDF save last L anyway
        out.Ls = [{L}, out.Ls(1:min(step, length(out.Ls)))];
    end

    out.Ks = [{nmout.K}, out.Ks];

    if isfield(opts.nm,'res_tol') && isnumeric(opts.nm.res_tol)
        out.res = [nmout.res(end), out.res];
    end

    if isfield(opts.nm,'rel_diff_tol') && isnumeric(opts.nm.rel_diff_tol)
        out.rc = [nmout.rc(end), out.rc];
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if eqn.LTV && not(opts.bdf.save_solution)
    out = rmfield(out, 'Ls');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
