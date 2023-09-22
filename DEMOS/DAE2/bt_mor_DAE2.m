function bt_mor_DAE2(problem, level, re, istest)
% Computes a standard ROM by implicitly solving the generalized Lyapunov
% equations for the equivalent projected system on the hidden manifold.
%
% Inputs:
% problem       either 'Stokes' or 'NSE' to choose the Stokes demo or the
%               linearized Navier-Stokes-Equation.
%               (optional, defaults to 'Stokes')
%
% level         discretization level 1 through 5
%               (optional, only used in 'NSE' case, default: 1)
%
% re            Reynolds number 300, 400, or 500
%               (optional, only used in 'NSE' case, default: 500)
%
% istest        flag to determine whether this demo runs as a CI test or
%               interactive demo
%               (optional, defaults to 0, i.e. interactive demo)
%
% Note that the 'NSE' option requires additional data available in a
% separate 270MB archive and at least the 5th discretization level needs a
% considerable amount of main memory installed in your machine.
%
% See:
% P. Benner, J. Saak, M. M. Uddin, Balancing based model reduction for
% structured index-2 unstable descriptor systems with application to flow
% control, Numerical Algebra, Control and Optimization 6 (1) (2016) 1–20.
% https://doi.org/10.3934/naco.2016.6.1

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% ADI tolerance and maximum iteration number
opts.adi.maxiter      = 350;
opts.adi.res_tol      = sqrt(eps);
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info         = 1;
opts.shifts.info      = 1;
opts.norm             = 'fro';

[oper, opts] = operatormanager(opts, 'dae_2');
%% Problem data
if nargin < 1
    problem = 'stokes';
end
if nargin < 2
    level   = 1;
end
if nargin < 3
    re      = 500;
end
if nargin < 4
    istest  = false;
end

problem = lower(problem);

switch problem
    case 'stokes'
        nin  = 5;
        nout = 5;
        nx   = 10;
        ny   = 10;
        [eqn.E_, eqn.A_, eqn.Borig, eqn.Corig] = ...
            stokes_ind2(nin, nout, nx, ny, opts);
        n = size(eqn.E_, 1);
        eqn.haveE = true;
        n_ode = trace(eqn.E_); % Stokes is FDM discretized, so so this is
        % the dimension of the velocity space
        eqn.manifold_dim = n_ode;
        eqn.B = eqn.Borig(1:n_ode, :);
        eqn.C = eqn.Corig(:, 1:n_ode);
    case 'nse'
        [eqn, K0primal, K0dual] = mess_get_NSE(re, level);
        n_ode = eqn.manifold_dim;
        n = size(eqn.E_, 1);
    otherwise
        mess_err(opts, ...
                 'illegal_input', ...
                 'input ''problem'' must be either ''NSE'' or ''Stokes''');
end
%%
eqn.type       = 'N';
% Activate stabilizing (Bernoulli) feedback
if strcmp(problem, 'nse')
    eqn.V      = -K0primal';
    eqn.U      = eqn.B;
    eqn.haveUV = true;
end

opts.shifts.num_desired = 6;
opts.shifts.num_Ritz    = 40;
opts.shifts.num_hRitz   = 40;
opts.shifts.method      = 'projection';

opts.shifts.b0 = ones(size(eqn.A_, 1), 1);

t_mess_lradi = tic;
outB = mess_lradi(eqn, opts, oper);
t_elapsed1 = toc(t_mess_lradi);
mess_fprintf(opts, ...
             'mess_lradi took %6.2f seconds \n', ...
             t_elapsed1);

if not(istest)
    figure('Name', 'Residual history (controllability)');
    semilogy(outB.res, 'LineWidth', 3);
    title('A X E^T + E X A^T = -BB^T');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
[mZ, nZ] = size(outB.Z);
mess_fprintf(opts, 'size outB.Z: %d x %d\n\n', mZ, nZ);

%%
eqn.type = 'T';
% Activate stabilizing (Bernoulli) feedback (for the dual system)
if strcmp(problem, 'nse')
    eqn.U      = -K0dual';
    eqn.V      = eqn.C';
    eqn.haveUV = true;
end

opts.shifts.num_desired = 6;
opts.shifts.num_Ritz    = 40;
opts.shifts.num_hRitz   = 40;
opts.shifts.method      = 'projection';

opts.shifts.b0 = ones(size(eqn.A_, 1), 1);

t_mess_lradi = tic;
outC = mess_lradi(eqn, opts, oper);
t_elapsed2 = toc(t_mess_lradi);
mess_fprintf(opts, ...
             'mess_lradi took %6.2f seconds \n', ...
             t_elapsed2);

if not(istest)
    figure('Name', 'Residual history (observability)');
    semilogy(outC.res, 'LineWidth', 3);
    title('A^T X E + E^T X A = -C^T C');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
[mZ, nZ] = size(outC.Z);
mess_fprintf(opts, 'size outC.Z: %d x %d\n\n', mZ, nZ);

%% Compute reduced system matrices
% Perform Square Root Method  (SRM)

% BT tolerance and maximum order for the ROM
t_SRM             = tic;
opts.srm.tol      = 1e-5;
opts.srm.max_ord  = 250;

% SRM verbosity
if istest
    opts.srm.info = 1;
else
    opts.srm.info = 2;
end

% The actual SRM
[TL, TR, hsv] = mess_square_root_method(eqn, opts, oper, outB.Z, outC.Z);

%%
ROM.A = TL' * (eqn.A_(1:n_ode, 1:n_ode) * TR);
ROM.B = TL' * eqn.B(1:n_ode, :);
ROM.C = eqn.C(:, 1:n_ode) * TR;

t_elapsed3 = toc(t_SRM);
mess_fprintf(opts, ...
             'computation of reduced system matrices took %6.2f seconds \n', ...
             t_elapsed3);

%%
t_eval_ROM = tic;
%% Evaluate the ROM quality
% while the Gramians are computed exploiting the DAE structure, due to the
% construction of the function handles we can not do so for the transfer
% function. Therefore we need to extend the matrices B and C and call the
% 'default' usfs for unstructured computation:
switch lower(problem)
    case 'stokes'
        eqn.B = eqn.Borig;
        eqn.C = eqn.Corig;
    case 'nse'
        n = size(eqn.A_, 1);
        eqn.B(n_ode + 1:n, :) = zeros(n - n_ode, size(eqn.B, 2));
        eqn.C(:, n_ode + 1:n) = zeros(size(eqn.C, 1), n - n_ode);
end
[oper, opts] = operatormanager(opts, 'default');

if istest
    opts.tf_plot.info = 0;
else
    opts.tf_plot.info = 2;
end

opts.tf_plot.fmin = -3;
opts.tf_plot.fmax =  4;

opts.tf_plot.type = 'sigma';

out = mess_tf_plot(eqn, opts, oper, ROM);
err = out.err;

t_elapsed4 = toc(t_eval_ROM);
mess_fprintf(opts, ...
             'evaluation of rom quality took %6.2f seconds \n', ...
             t_elapsed4);

%%
if istest
    if max(err) >= opts.srm.tol
        mess_err(opts, 'TEST:accuracy', 'unexpectedly inaccurate result');
    end
else
    figure('Name', 'Computed Hankel singular values');
    semilogy(hsv, 'LineWidth', 3);
    title('Computed Hankel singular values');
    xlabel('index');
    ylabel('magnitude');
end
%%
mess_fprintf(opts, ...
             ['\nComputing open loop step response of original and reduced ' ...
              'order systems and time domain MOR errors\n']);
open_step(eqn, ROM.A, ROM.B, ROM.C, problem, istest);
%%
mess_fprintf(opts, '\nComputing ROM based feedback\n');
if exist('care', 'file')
    [~, ~, Kr] = care(ROM.A, ROM.B, ROM.C' * ROM.C, eye(size(ROM.B, 2)));
else
    Y = care_nwt_fac([], ROM.A, ROM.B, ROM.C, 1e-12, 50);
    Kr = (Y * ROM.B)' * Y;
end
K = [Kr * TL' * eqn.E_(1:n_ode, 1:n_ode), zeros(size(Kr, 1), n - n_ode)];
%%
mess_fprintf(opts, ...
             ['\nComputing closed loop step response of original and ' ...
              'reduced order systems and time domain MOR errors\n']);
closed_step(eqn, ROM.A, ROM.B, ROM.C, problem, K, Kr, istest);
