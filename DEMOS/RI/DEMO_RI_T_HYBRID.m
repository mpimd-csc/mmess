function DEMO_RI_T_HYBRID(istest)
% Computes the solution of the Hinf Riccati equation for a random generated
% generalized system. The computations are done using the Newton method for
% the LQG step and afterwards the RADI. Afterwards, the real residual
% norm is shown and compared to the set tolerance.
%
% Input:
% istest  decides whether the function runs as an interactive demo or a
%         continuous integration test. (optional; defaults to 0, i.e.
%         interactive demo)
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
if nargin < 1
    istest = false;
end

%% Construction of system data.
if exist('OCTAVE_VERSION', 'builtin')
    rand('seed', 1.0); %#ok<RAND>
    eqn.A_ = rand(500) - 250 * eye(500);

    rand('seed', 2.0); %#ok<RAND>
    eqn.B2 = rand(500, 2);
    rand('seed', 3.0); %#ok<RAND>
    B1 = rand(500, 2);

    rand('seed', 4.0); %#ok<RAND>
    eqn.C2 = rand(3, 500);
    rand('seed', 5.0); %#ok<RAND>
    C1 = rand(3, 500);
else
    rng(1.0);
    eqn.A_ = rand(500) - 250 * eye(500);

    rng(2.0);
    eqn.B2 = rand(500, 2);
    rng(3.0);
    B1 = rand(500, 2);

    rng(4.0);
    eqn.C2 = rand(3, 500);
    rng(5.0);
    C1 = rand(3, 500);
end

gam = 5; % Scaling term for disturbances.

%% Set operator.
opts = struct();
[oper, opts] = operatormanager(opts, 'default');

%% Construction of options struct.
% ADI settings.
opts.adi.maxiter          = 200;
opts.adi.res_tol          = 1.0e-12;
opts.adi.rel_diff_tol     = 0;
opts.adi.info             = 1;
opts.adi.compute_sol_fac  = true;
opts.adi.accumulateK      = false;
opts.adi.accumulateDeltaK = false;

% Shift options.
opts.shifts.num_desired   = 5;
opts.shifts.method = 'projection';

% % NM settings.
opts.nm.maxiter       = 50;
opts.nm.res_tol       = 1.0e-10;
opts.nm.rel_diff_tol  = 1.0e-12;
opts.nm.info          = 1;
opts.nm.linesearch    = false;
opts.nm.accumulateRes = true;

% NM projection settings.
opts.nm.projection      = [];
opts.nm.projection.freq = 0;
opts.nm.res.maxiter     = 10;
opts.nm.res.tol         = 1.0e-06;
opts.nm.res.orth        = true;

% RADI settings.
opts.radi.maxiter       = opts.adi.maxiter;
opts.radi.res_tol       = opts.nm.res_tol;
opts.radi.rel_diff_tol  = 1.0e-16;
opts.radi.info          = 1;

% RI settings.
opts.ri.riccati_solver = 'radi';
opts.ri.lqg_solver     = 'newton';
opts.ri.maxiter        = 10;
opts.ri.res_tol        = 1.0e-10;
opts.ri.rel_diff_tol   = 1.0e-16;
opts.ri.compres_tol    = 1.0e-16;
opts.ri.info           = 1;

% global options
opts.norm           = 2;

% %% Call Riccati iteration with Newton solver.
t_RI_call = tic;
eqn.type = 'T';
eqn.B1   = 1 / gam * B1;
eqn.C1   = C1;
[out, eqn, opts, ~] = mess_lrri(eqn, opts, oper);
t_elapsed = toc(t_RI_call);
mess_fprintf(opts, 'mess_lrri took %6.2f seconds \n', t_elapsed);

%% Compute real residuals.
abserr = norm(eqn.A_' * (out.Z * out.Z') + (out.Z * out.Z') * eqn.A_ + ...
              (out.Z * out.Z') * (1 / gam^2 * (B1 * B1') - ...
                                  eqn.B2 * eqn.B2') * ...
              (out.Z * out.Z') + C1' * C1, 2);
relerr = abserr / norm(C1 * C1', 2);
mess_fprintf(opts, '\nset tolerance vs. real residual: %e | %e\n', ...
             opts.ri.res_tol, relerr);

if istest
    mess_assert(opts, relerr < opts.ri.res_tol, ...
                'TEST:accuracy', 'unexpectedly inaccurate result');
end
