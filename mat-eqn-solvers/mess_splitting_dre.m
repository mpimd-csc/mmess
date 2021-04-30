 function [out, eqn, opts, oper] = mess_splitting_dre(eqn, opts, oper)
%% function [out, eqn, opts, oper] = mess_splitting_dre(eqn, opts, oper)
%   LDL^T-factored splitting schemes for diffenential Riccati equations
%   E*d/dt X(t)*E' = -B*B' - E*X(t)*A' - A*X(t)*E' + E*X(t)*C'*C*X(t)*E' (N)
%   E'*d/dt X(t)*E = -C'*C - E'*X(t)*A - A'*X(t)*E + E'*X(t)*B*B'*X(t)*E (T)
%   backward in time.
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
%   out     struct contains solutions for every time step
%
% Input fields in struct eqn:
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   eqn.L0      dense (n x m4) matrix, L in initial LDL^T factorization
%
%   eqn.D0      dense (m4 x m4) matrix, D in initial LDL^T factorization
%
%   eqn.type    possible values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional, default 'N')
%
%   eqn.haveE   possible values: 0, 1, false, true
%               if haveE = 0: matrix E is assumed to be the identity
%               (optional)
%
%   eqn.LTV     possible  values: 0, 1, false, true
%               indicates autonmous (false) or
%               non-autonomous (true) differential
%               Riccati equation
%               (optional, default: 0)
%
%   eqn.Rinv    dense (m1 x m1) matrix implementing the term
%               E'*X(t)*B*Rinv*B'*X(t)*E (T)
%               (Rinv is not used in computation of out.Ks)
%
%   non-autonomous (LTV), NOTE: A(t)*A(s)=A(s)*A(t) is required. This is
%   not checked in the code and has to be ensured by the user.
%
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices. For the second order types these are given
%   implicitly by the M, D, K matrices stored in eqn.M_, eqn.E_ and eqn.K_,
%   respectively.
%
% Input fields in struct opts:
%   opts.splitting.time_steps   possible values: (N x 1) array
%                               array containing the equidistant time steps
%                               If splitting.adaptive is set, it should
%                               specify the time interval [t0, tend]
%
%   opts.splitting.order        possible values: 1, 2, 3, 4, ...
%                               use splitting scheme of this order
%                               1 is Lie splitting, 2 is Strang splitting
%                               (default), 3- requires splitting.additive
%                               to be true
%                               (optional, default 2)
%
%   opts.splitting.additive     possible values: 0, 1, false, true
%                               use (1) or don't use (0) the additive
%                               schemes. If splitting.order > 2 this is set
%                               to 1 automatically
%                               (optional)
%
%   opts.splitting.symmetric    possible values: 0, 1, false, true
%                               use symmetric (1) or unsymmetric (0)
%                               additive schemes
%                               (optional, default 0)
%
%   opts.splitting.trunc_tol    possible values: scalar > 0
%                               column compression tolerance
%                               (optional, default eps * size(A, 1))
%
%   opts.splitting.trunc_info   possible values: 0, 1, false, true
%                               verbose mode for column compression
%                               (optional, default: 0)
%
%   opts.splitting.quadrature.type
%                               possible values: 'gauss', 'clenshawcurtis',
%                                                'equidistant', 'adaptive'
%                               type of quadrature to use for approximation
%                               of the integral term
%                               (optional, default 'gauss')
%
%   opts.splitting.quadrature.order
%                               possible values: positive integer
%                               order of quadrature  used for approximation
%                               of the integral term
%                               (optional, default based on method)
%
%   opts.splitting.quadrature.tol
%                               possible values: positive real value
%                               tolerance to use for adaptive quadrature
%                               (optional, default 1e-4)
%
%   opts.splitting.intermediates
%                               possible values: 0, 1, false, true
%                               store intermediate approximations (1)
%                               or only return final L, D (0)
%                               (optional, default 1)
%
%   opts.splitting.info         possible values: 0, 1, 2, false, true
%                               turn on (1) or off (0) the status output at
%                               every time step. 2 prints only the time
%                               steps.
%                               (optional, default 2)
%
%
%  The structure opts.exp_action contains parameters for how to
%  compute matrix exponential actions on block matrices; see
%  private/mess_exp_action.m for specifications
%
%
%
%
%
% % Adaptivity, NOTE: not implemented, future functionality
%
%   opts.splitting.adaptive     possible values: 0, 1, false, true, struct
%                               use adaptive time stepping (1, struct) or
%                               not (0)
%                               (optional, default 0)
%
%   opts.splitting.adaptive.controller  possible values: 'PI',0,'deadbeat'
%                                       the type of adaptive controller; PI
%                                       (with unoptimized parameters) or
%                                       deadbeat (0, 'deadbeat')
%                                       (optional, default 'PI')
%
%   opts.splitting.adaptive.initial_timestep possible values: scalar > 0
%                                            guess for good initial time
%                                            step
%                                            (optional)
%
%   opts.splitting.adaptive.epus        possible values: 0, 1, false, true
%                                       use the error per unit step (1)
%                                       strategy or error per step (0)
%                                       (optional, default 1)
%
%   opts.splitting.adaptive.reuse_IQ    possible values: 0, 1, false, true
%                                       reuse most data for computing I_Q
%                                       after the step size changes (1) or
%                                       not (0). Usually saves a lot of
%                                       computation time, but sometimes
%                                       leads to many convergence failures
%                                       (optional, default false)
%
% Output fields in struct out:
%
%   out.Ls      cell array with solution factor L for every time step
%               if opts.splitting.intermediates = 1, otherwise the final
%               factor
%
%   out.Ds      cell array with solution factor D for every time step
%               if opts.splitting.intermediates = 1, otherwise the final
%               factor
%
%   out.Ks      cell array with feedback term K for every time step
%
%   out.ranks   Nt x 1 array with numerical ranks of approximations for
%               every time step (saved also if opts.splitting.intermediates
%               = 0)
%
%   out.time    the total computation time (wall clock)
%
%   out.timeIQ  the computation time (wall clock) required for the I_Q
%               integral approximation
%
%
% If optional input arguments are missing they may be set to default values
% and a 'MESS:control_data' warning is printed. to turn warnings off use
% warning('OFF', 'MESS:control_data')
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for splitting  control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isfield(opts,'splitting')) || not(isstruct(opts.splitting))
    error('MESS:control_data', ['splitting control structure opts.' ...
                                'splitting missing.']);
end % Single fields are checked below or inside subfunctions

if not(isfield(opts.splitting,'time_steps')) ...
    && (not(isfield(opts.splitting,'adaptive')) || not(opts.splitting.adaptive))
    error('MESS:control_data', ['opts.splitting.time_steps is missing,' ...
                    'and adaptive time stepping has not been requested.']);
end
opts.t0 = opts.splitting.time_steps(1);

if not(isfield(opts.splitting, 'order')),  opts.splitting.order = 2; end

if rem(opts.splitting.order, 1) || opts.splitting.order < 1
    error('MESS:control_data', ['opts.splitting.order has an invalid ', ...
                                'value.']);
end


if not(isfield(opts.splitting, 'additive')),  opts.splitting.additive = 0; end
if opts.splitting.order  >  2, opts.splitting.additive = 1; end

if not(isfield(opts.splitting, 'symmetric')) && opts.splitting.additive == 1
    opts.splitting.symmetric = 0;
end

if opts.splitting.additive && opts.splitting.symmetric ...
                           && rem(opts.splitting.order, 2)
    error('MESS:control_data', ['opts.splitting.order must be a ' ...
                'multiple of 2 to use a symmetric scheme.']);
end

if not(isfield(opts.splitting, 'quadrature'))
    error('MESS:control_data', ['Need to specify ' ...
          'opts.splitting.quadrature struct.']);
end

if not(isfield(opts.splitting.quadrature, 'type'))
        warning('MESS:control_data',['Unspecified quadrature type. ' ...
                'Setting opts.splitting.quadrature.type=''gauss''']);
    opts.splitting.quadrature.type = 'gauss';
elseif not(ismember(opts.splitting.quadrature.type, ...
                {'gauss', 'clenshawcurtis', 'equidistant', 'adaptive'}))
    warning('MESS:control_data',['The specified quadrature type is ' ...
                        'not supported.\n setting ' ...
                        'opts.splitting.quadrature.type=''gauss''']);
    opts.splitting.quadrature.type = 'gauss';
end

if not(isfield(opts.splitting.quadrature, 'order'))
    s = opts.splitting.order;
    if s == 1 % Lie
        opts.splitting.quadrature.order = 3;
    end
    if s == 2 % Strang
        opts.splitting.quadrature.order = 5;
    end
    if opts.splitting.additive
        if opts.splitting.symmetric
            opts.splitting.quadrature.order = 2*s + 3;
        else
            if mod(s, 2) == 0 % even
                opts.splitting.quadrature.order = s + 3;
            else % odd
                opts.splitting.quadrature.order = s + 2;
            end
        end
    end

end

if strcmp(opts.splitting.quadrature.type, 'clenshawcurtis')
    opts.splitting.quadrature.order = opts.splitting.quadrature.order ...
        + rem(opts.splitting.quadrature.order, 2);
end

if strcmp(opts.splitting.quadrature.type, 'adaptive') && ...
   not(isfield(opts.splitting.quadrature, 'tol'))
        warning('MESS:control_data',['Using adaptive quadrature, but ' ...
                        'tolerance unspecified. Setting ' ...
                        'opts.splitting.quadrature.tol=1e-4']);
	opts.splitting.quadrature.tol=1e-4;
end

if not(isfield(opts.splitting, 'intermediates'))
    opts.splitting.intermediates = 1;
end


if not(isfield(opts.splitting,'info')),  opts.splitting.info = 2; end

if isfield(opts,'LDL_T') && opts.LDL_T==0
    warning('MESS:control_data',['The splitting code only supports ' ...
                        'LDL_T solutions.\n Setting opts.LDL_T=1']);
end
opts.LDL_T = 1;





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for time step adaptiveness control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isfield(opts.splitting, 'adaptive'))  ...
    || not(isstruct(opts.splitting.adaptive))
    opts.splitting.adaptive = 0;
end

if opts.splitting.adaptive ~= 0

    error('MESS:missing_feature', ['Time step adaptivity has not yet been'
                                    'implemented.']);

    if not(isfield(opts.splitting.adaptive, 'controller')) %#ok<UNRCH>
        opts.splitting.adaptive.controller = 'PI';
    end

    if not(isfield(opts.splitting.adaptive, 'initial_timestep'))
        tend = opts.splitting.time_steps(end);
        % Very crude and probably overly pessimistic guess
        opts.splitting.adaptive.initial_timestep = (tend - t0) / 1000;
    end

    if not(isfield(opts.splitting.adaptive, 'epus'))
        opts.splitting.adaptive.epus = 1;
    end

    if not(isfield(opts.splitting.adaptive, 'reuse_IQ'))
        opts.splitting.adaptive.reuse_IQ = 0;
    end

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn,'type'))
    eqn.type='N';
    warning('MESS:control_data',['Unable to determine type of equation.'...
        'Falling back to type ''N''']);
elseif (eqn.type~='N') && (eqn.type~='T')
    error('MESS:equation_type', ['Equation type must be either ''T'' ', ...
                                 'or ''N''']);
end

if not(isfield(eqn, 'LTV'))
    eqn.LTV = 0;
end

if eqn.LTV % Instantiate matrices at first time step
    if isfield(oper, 'eval_matrix_functions')
        [eqn, opts, oper] = oper.eval_matrix_functions(eqn, opts, oper, ...
                                                opts.splitting.time_steps(1));
    else
        error('MESS:missing_feature', ['The function eval_matrix_functions is ', ...
                                       'required for LTV problems, but it has ', ...
                                       'not yet been implemented for this set ', ...
                                       'of USFS functions']);
    end
end

[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');

if not(result)
    error('MESS:control_data', ['system data is not completely ', ...
                                'defined or corrupted']);
end

% If the user sends in sparse B and C, make them dense
if not(eqn.LTV) && issparse(eqn.B)
    eqn.B = full(eqn.B);
    eqn.C = full(eqn.C);
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
    eqn.L0 = zeros(oper.size(eqn, opts), 1);
end

if not(isfield(eqn, 'D0')) || not(isnumeric(eqn.D0))
    warning('MESS:control_data', ...
            ['Initial condition factor D0 is not defined or corrupted.',...
             ' Setting it to the identity matrix.']);
    eqn.D0 = eye(size(eqn.L0, 2));
end


% Extra check of general splitting data, which has to be run after the
% operator structure has been initialized. It cannot be initialized before
% the main splitting checks, because in the LTV case the initialization
% needs the first time step, which exists in opts.splitting.time_steps...
if not(isfield(opts.splitting,'trunc_tol'))
    opts.splitting.trunc_tol = eps * oper.size(eqn, opts);
end
if not(isfield(opts.splitting,'trunc_info'))
    opts.splitting.trunc_info = 0;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All checks done. Here comes the real work!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We solve E'*d/dt X(t)*E = C'*C + E'*X(t)*A + A'*X(t)*E
%                           - E'*X(t)*B*B'*X(t)*E (eqn.type == 'T')
% forward in time instead and then flip everything.
t = opts.splitting.time_steps;
% In the LTV case, this means we have to evaluate the matrix functions at
% tend - t rather than t.
if eqn.LTV
    opts.splitting.eval_matrix_functions = @(eqn, opts, oper, s) ...
        oper.eval_matrix_functions(eqn, opts, oper, t(end) - s);
else
    opts.splitting.eval_matrix_functions = @(eqn, opts, oper, s) ...
        mess_do_nothing(eqn, opts, oper);
end


Nt = length(t) - 1;
h = (t(end) - t(1)) / Nt;

L0 = eqn.L0;
D0 = eqn.D0;

ms = zeros(Nt, 1);
ms(1) = size(L0, 2);

if opts.splitting.intermediates
    Ls{Nt+1} = [];
    Ls{1} = L0;
    Ds{Nt+1} = [];
    Ds{1} = D0;
else
    Ls = L0;
    Ds = D0;
end

Ks{Nt+1} = [];
[eqn, opts, oper] = ...
    opts.splitting.eval_matrix_functions(eqn, opts, oper, t(1));
if eqn.type == 'T'
    Ks{1} = (oper.mul_E(eqn, opts, 'T', L0*(D0*(L0'*eqn.B)), 'N'))';
elseif eqn.type == 'N'
    Ks{1} = (oper.mul_E(eqn, opts, 'T', L0*(D0*(L0'*eqn.C')), 'N'))';
end
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up method coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = opts.splitting.order;
if s == 1 % Lie
    as = 1;
end
if s == 2 % Strang
    as = 1/2;
end
if opts.splitting.additive
    if opts.splitting.symmetric
        switch s/2
            case 1
                gamma = 1/2;
            case 2
                gamma = [-1/6, 2/3];
            case 3
                gamma = [1/48; -8/15; 81/80];
            case 4
                gamma = [-1/720; 8/45; -729/560; 512/315];
            otherwise
                gamma = compute_additive_coefficients(s, true);
        end
        as = 1./(1:s/2);
    else
        switch s
            case 1
                gamma = 1;
            case 2
                gamma = [-1, 2];
            case 3
                gamma = [1/2, -4, 9/2];
            case 4
                gamma = [-1/6, 4, -27/2, 32/3];
            otherwise
                gamma = compute_additive_coefficients(s, false);
        end
        as = 1./(1:s);
    end
end

% If autonomous, we may approximate the integral once and for all
if not(eqn.LTV)
    [IQL, IQD] = IQ(eqn, opts, oper, 0, h, as);
%     % Extra output for time-step adaptivity (not yet implemented):
%     out.IQLas = IQLas;
%     out.IQDas = IQDas;
%     out.Las = Las;
%     out.xj = xj;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1: Nt
    if opts.splitting.info > 0
        fprintf('Step: %d of %d, time = %d \n', j, Nt, t(j));
    end
    if opts.splitting.intermediates
        Lold = Ls{j}; Dold = Ds{j};
    else
        Lold = Ls; Dold = Ds;
    end

    if s == 1 % Lie
        [Lnew, Dnew] = expG(eqn, opts, oper, h, Lold, Dold, t(j));
        if eqn.LTV
            [IQL, IQD] = IQ(eqn, opts, oper, t(j), h, as);
        end
        [Lnew, Dnew, eqn, opts, oper] ...
        = expF(eqn, opts, oper, h, IQL{1}{1}, IQD{1}{1}, Lnew, Dnew, t(j));
    elseif s == 2 && not(opts.splitting.additive) % Strang
        if eqn.LTV
            [IQL, IQD] = IQ(eqn, opts, oper, [t(j), t(j) + h/2], h, as);
        end
        [Lnew, Dnew, eqn, opts, oper] ...
            = expF(eqn, opts, oper, ...
                    h/2, IQL{1}{1}, IQD{1}{1}, Lold, Dold, t(j));
        [Lnew, Dnew] = expG(eqn, opts, oper, h, Lnew, Dnew, t(j));
        if eqn.LTV
            [Lnew, Dnew, eqn, opts, oper] ...
                = expF(eqn, opts, oper, ...
                       h/2, IQL{2}{1}, IQD{2}{1}, Lnew, Dnew, t(j) + h/2);
        else
            [Lnew, Dnew, eqn, opts, oper] ...
            = expF(eqn, opts, oper, h/2, IQL{1}{1}, IQD{1}{1}, Lnew, Dnew);
        end
    elseif opts.splitting.additive
        if opts.splitting.symmetric, nmax = s/2; else, nmax = s; end

        tj = t(j);

        if opts.splitting.symmetric
            if not(eqn.LTV)
                IQL1 = IQL{1};
                IQD1 = IQD{1};

                parfor n = 1:nmax
                    % Local copy of opts for parallelization
                    opts_par = opts;
                    L1{n} = Lold; D1{n} = Dold;
                    L2{n} = Lold; D2{n} = Dold;

                    for k = 1:n
                        [L2{n}, D2{n}, ~, opts_par, ~] ...
                            = expF(eqn, opts_par, oper, ...
                                   h / n, IQL1{n}, IQD1{n}, L2{n}, D2{n});
                        [L2{n}, D2{n}] ...
                            = expG(eqn, opts, oper, h / n, L2{n}, D2{n});
                        [L1{n}, D1{n}] ...
                            = expG(eqn, opts, oper, h / n, L1{n}, D1{n});
                        [L1{n}, D1{n}, ~, opts_par, ~] ...
                            = expF(eqn, opts_par, oper, ...
                                   h / n, IQL1{n}, IQD1{n}, L1{n}, D1{n});
                    end
                end
            else % Time-varying
                parfor n = 1:nmax
                    % Local copy of opts for parallelization
                    opts_par = opts;
                    L1{n} = Lold; D1{n} = Dold;
                    L2{n} = Lold; D2{n} = Dold;

                    % New name, because otherwise Matlab complains that the
                    % temporary variable is used after the parfor-loop.
                    % (Even though it isn't.)
                    [IQLpar, IQDpar] = ...
                        IQ(eqn, opts_par, oper, tj + h*(0:n-1)/n, h, 1/n);

                    for k = 1:n
                        tt = tj + (k-1) / n * h;
                        [L2{n}, D2{n}, ~, opts_par, ~] ...
                            = expF(eqn, opts_par, oper, ...
                                   h / n, IQLpar{k}{1}, IQDpar{k}{1}, ...
                                   L2{n}, D2{n}, tt);
                        [L2{n}, D2{n}] = expG(eqn, opts, oper, h / n, ...
                                              L2{n}, D2{n}, tt);
                        [L1{n}, D1{n}] = expG(eqn, opts, oper, h / n, ...
                                              L1{n}, D1{n}, tt);
                        [L1{n}, D1{n}, ~, opts_par, ~] ...
                            = expF(eqn, opts_par, oper, ...
                                   h / n, IQLpar{k}{1}, IQDpar{k}{1}, ...
                                   L1{n}, D1{n}, tt);
                    end
                end
            end

            % Stack the results
            Lnew = [cell2mat(L1), cell2mat(L2)];

            for n = 1:nmax
                D1{n} = gamma(n)*D1{n};
                D2{n} = gamma(n)*D2{n};
            end
            Dnew = blkdiag(D1{:}, D2{:});

            [Lnew, Dnew] = mess_column_compression(Lnew, 'N', Dnew, ...
                opts.splitting.trunc_tol, opts.splitting.trunc_info);

        else % Asymmetric case
            if not(eqn.LTV)
                IQL1 = IQL{1};
                IQD1 = IQD{1};

                parfor n = 1:nmax
                    % Local copy of opts for parallelization
                    opts_par = opts;
                    L1{n} = Lold; D1{n} = Dold;

                    for k = 1:n
                        [L1{n}, D1{n}] ...
                            = expG(eqn, opts, oper, h / n, L1{n}, D1{n});
                        [L1{n}, D1{n}, ~, opts_par, ~] ...
                            = expF(eqn, opts_par, oper, ...
                                   h / n, IQL1{n}, IQD1{n}, L1{n}, D1{n});
                    end
                end
            else % Time-varying case
                parfor n = 1:nmax
                    % Local copy of opts for parallelization
                    opts_par = opts;
                    L1{n} = Lold; D1{n} = Dold;

                    % New name, because otherwise Matlab complains that the
                    % temporary variable is used after the parfor-loop.
                    % (Even though it isn't.)
                    [IQLpar, IQDpar] = ...
                        IQ(eqn, opts_par, oper, tj + h*(0:n-1)/n, h, 1/n);

                    for k = 1:n
                        tt = tj + (k-1) / n * h;
                        [L1{n}, D1{n}] = expG(eqn, opts, oper, h / n, ...
                                              L1{n}, D1{n}, tt);
                        [L1{n}, D1{n}, ~, opts_par, ~] ...
                            = expF(eqn, opts_par, oper, ...
                                   h / n, IQLpar{k}{1}, IQDpar{k}{1}, ...
                                   L1{n}, D1{n}, tt);
                    end
                end
            end

            % Stack the results
            Lnew = cell2mat(L1);

            for n = 1:nmax
                D1{n} = gamma(n)*D1{n};
            end
            Dnew = blkdiag(D1{:});

            [Lnew, Dnew] = mess_column_compression(Lnew, 'N', Dnew, ...
                opts.splitting.trunc_tol, opts.splitting.trunc_info);
        end
    end

    if opts.splitting.intermediates
        Ls{j+1} = Lnew; Ds{j+1} = Dnew;
    else
        Ls = Lnew; Ds = Dnew;
    end
    ms(j+1) = size(Dnew, 1);

    % Store feedback term K as well. We have
    %  (T): K = B' LDL' E = (E'*LDL'*B)' and
    %  (N): K = C LDL' E = (E'*LDL'*C')'
    % The extra transpose is because we can only multiply with E from the
    % left.
    if eqn.type == 'T'
        Ks{j+1} = ...
          (oper.mul_E(eqn, opts, 'T', Lnew*(Dnew*(Lnew'*eqn.B)), 'N'))';
    elseif eqn.type == 'N'
        Ks{j+1} = ...
          (oper.mul_E(eqn, opts, 'T', Lnew*(Dnew*(Lnew'*eqn.C')), 'N'))';
    end

end

% Flip to go backward rather than forward in time
out.Ls = fliplr(Ls);
out.Ds = fliplr(Ds);
out.ms = flipud(ms);
out.Ks = fliplr(Ks);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);


 end