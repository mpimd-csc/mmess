function out = LQR_LTV_smallscale_splitting(method, istest)

% Computes the optimal feedback via low-rank splitting schemes [1, 2] for
% a small-scale LTV system, where the matrices are time-varying in the
% sense that A(t) = a(t) \hat{A}.
%
% Inputs:
%
% method      choice of splitting scheme; structure with fields 'order',
%             'additive' and 'symmetric'
%             (optional, defaults to order = 2, additive = false, symmetric
%             = false)
%             NOTE: the additive schemes typically need a running parallel
%             pool in order to be competitive
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
%
% References:
%
% [1] T. Stillfjord, Low-rank second-order splitting of large-scale
%     differential Riccati equations, IEEE Trans. Autom. Control, 60
%     (2015), pp. 2791-2796. https://doi.org/10.1109/TAC.2015.2398889
%
% [2] T. Stillfjord, Adaptive high-order splitting schemes for large-scale
%     differential Riccati equations, Numer. Algorithms,  (2017).
%     https://doi.org/10.1007/s11075-017-0416-8

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin < 1
    method.order = 2;
    method.additive = false;
    method.symmetric = false;
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
e = ones(Nx, 1);
A_1D = 1 / dx^2 * spdiags([e, -2 * e, e], -1:1, Nx, Nx);

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

% R^{-1}, inverse of weighting factor for input in cost functional
eqn.Rinv = 1;
opts = struct();
[oper, opts] = operatormanager(opts, 'default');

eqn.LTV = true; % Specify that this is a time-varying problem

% Time interval [0, 0.1] and 100 time steps
t0 = 0;
tend = 0.1;
Nt = 100;

%% General splitting parameters
opts.splitting.time_steps = linspace(t0, tend, Nt + 1);
opts.splitting.order = method.order;
opts.splitting.additive = method.additive;
opts.splitting.symmetric = method.symmetric;
opts.splitting.info = 2;
opts.splitting.intermediates = true;

opts.splitting.trunc_tol = eps;

% Quadrature (for integral terms) parameters
opts.splitting.quadrature.type = 'adaptive';
opts.splitting.quadrature.tol = 1e-8;

%% Matrix exponential actions
opts.exp_action.method = 'LTV';
opts.exp_action.tol = 1e-8;

%% Compute the approximation
t_mess_splitting_dre = tic;
[out, ~, opts, ~] = mess_splitting_dre(eqn, opts, oper);
t_elapsed = toc(t_mess_splitting_dre);
mess_fprintf(opts, 'mess_splitting_dre took %6.2f seconds \n', t_elapsed);

%%
if not(istest)
    t = opts.splitting.time_steps;
    figure;
    plot(t, out.ms, 'LineWidth', 3);
    title('Ranks of approximations over time');

    y = zeros(1, length(out.Ks));
    for i = 1:length(out.Ks)
        y(i) = out.Ks{i}(1, 1);
    end

    figure;
    plot(t, y, 'LineWidth', 3);
    title('evolution of component (1,1) of the optimal feedback');
else

    if abs(norm(out.Ds{1}) / 6.078945091766749e-05 - 1) >= 1e-10
        mess_err(opts, 'TEST:accuracy', 'unexpectedly inaccurate result');
    end
end
