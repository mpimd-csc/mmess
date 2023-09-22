function out = LQR_LTV_smallscale_BDF(k, istest)

% Computes the optimal feedback via a BDF method for a small-scale LTV
% system, where the matrices are time-varying in the sense that
% A(t) = a(t) \hat{A}.
%
% Inputs:
%
% k           k-step BDF method
%             possible values: 1, ..., 4
%             (optional, defaults to 2)
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin < 1
    k = 2;
end
if nargin < 2
    istest = false;
end

Nx = 10;
Nx2 = Nx^2;

% Finite-difference approximation of Laplacian on [0,1]^2 with Dirichlet
% boundary conditions
xx = linspace(0, 1, Nx + 2);
x = xx(2:end - 1);
dx = 1 / (Nx + 1);
ev = ones(Nx, 1);
A_1D = 1 / dx^2 * spdiags([ev, -2 * ev, ev], -1:1, Nx, Nx);

[X, Y] = meshgrid(x, x);
I_1D = eye(Nx, Nx);
A_2D = kron(A_1D, I_1D) + kron(I_1D, A_1D);
M = speye(Nx2, Nx2); % mass matrix = I in this case

% 3 inputs, control on the squares [j/4, j/4 + 1/8] x [j/4, j/4 + 1/8]
Nb = 3;
B = zeros(Nx2, Nb);
for j = 1:Nb
    B(:, j) = reshape((X > j / 4) & (X < j / 4 + 1 / 8) & ...
                      (Y > j / 4) & (Y < j / 4 + 1 / 8), ...
                      Nx2, 1);
end

% 1 output, average of all states
C = 1 / Nx2 * ones(1, Nx2);

% Time dependency through factors
alpha = @(t) 1 + 10 * sin(2 * pi * t); % A
mu = @(t) 2 + 7.1 * sin(2 * pi * t);   % M
dmu = @(t) 7.1 * 2 * pi * cos(2 * pi * t); % dM/dt
Beta = @(t) 3 + cos(t);          % B
Gamma = @(t) 1 - min(t, 1);      % C

% Alternative, autonomous case:
% alpha = @(t) 1;        % A
% mu = @(t) 1;           % M
% dmu = @(t) 0;          % dM/dt
% beta = @(t) 1;         % B
% gamma = @(t) 1;        % C

At = @(t) alpha(t) * A_2D;
Mt = @(t) mu(t) * M;
dMt = @(t) dmu(t) * M;
Bt = @(t) Beta(t) * B;
Ct = @(t) Gamma(t) * C;

eqn.A_time = At;
eqn.E_time = Mt;
eqn.dt_E_time = dMt;
eqn.B_time = Bt;
eqn.C_time = Ct;

eqn.haveE = true;

eqn.type = 'T';

% LDLT-factorization of initial condition P(0) = 0
L0 = zeros(Nx2, 1);
eqn.L0 = L0;
eqn.D0 = eye(size(L0, 2));

eqn.LTV = true; % Specify that this is a time-varying problem

% Time interval [0, 0.1] and 100 time steps
t0 = 0;
tend = 0.1;
Nt = 1000;

% Set up and initialize operator
opts = struct;
[oper, opts] = operatormanager(opts, 'default');

[eqn, opts, oper] = oper.eval_matrix_functions(eqn, opts, oper, tend);

%% General BDF parameters

%%
opts.norm = 'fro';

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.res_tol = 1e-14;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;
opts.adi.compute_sol_fac = true;
opts.cc_info = 0;

%%
% Heuristic shift parameters via basic Arnoldi
opts.shifts.num_desired = 7;
opts.shifts.num_Ritz = 50;
opts.shifts.num_hRitz = 25;
opts.shifts.method = 'heur';

% opts.shifts.b0=ones(n,1);

%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 8;
opts.nm.res_tol = 1e-10;
opts.nm.rel_diff_tol = 1e-16;
opts.nm.info = 0;
opts.norm = 'fro';
opts.nm.accumulateRes = true;
opts.nm.linesearch = true;

%%
% BDF parameters
opts.bdf.time_steps = linspace(t0, tend, Nt + 1);
opts.bdf.step = k;
opts.bdf.info = 1;
opts.bdf.save_solution = 1;
opts.bdf.startup_iter = 7;

%% Compute the approximation
t_mess_bdf_dre = tic;
[out, ~, opts, ~] = mess_bdf_dre(eqn, opts, oper);
t_elapsed = toc(t_mess_bdf_dre);
mess_fprintf(opts, 'mess_bdf_dre took %6.2f seconds \n', t_elapsed);

%%
if not(istest)
    t = opts.bdf.time_steps;

    y = zeros(1, length(out.Ks));
    for i = 1:length(out.Ks)
        y(i) = out.Ks{i}(1, 1);
    end

    figure;
    plot(t, y, 'LineWidth', 3);
    title('evolution of component (1,1) of the optimal feedback');
else

    if abs(norm(out.Ds{1}) / 6.079051242083189e-05 - 1) >= 1e-10
        mess_err(opts, 'TEST:accuracy', 'unexpectedly inaccurate result');
    end
end
