function DEMO_RI_GE_DAE2(istest)
% Computes the solution of the Hinf Riccati equation for a random generated
% DAE index-2 system. The computations are done first time with the Newton
% method and the second time with the RADI. Afterwards, the explicitly
% projected equations is used to verify the results.
%
% Input:
% istest  decides whether the function runs as an interactive demo or a
%         continuous integration test. (optional; defaults to 0, i.e.
%         interactive demo)
%

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
if nargin<1, istest=0; end

%% Construction of system data.
if exist('OCTAVE_VERSION', 'builtin')
    rand('seed', 1.0); %#ok<RAND>
    A = rand(500) - 250 * eye(500);
    rand('seed', 2.0); %#ok<RAND>
    M = rand(500);
    M = M' * M;
    rand('seed', 3.0); %#ok<RAND>
    J = rand(100, 500);

    rand('seed', 4.0); %#ok<RAND>
    eqn.B2 = rand(500, 2);
    rand('seed', 5.0); %#ok<RAND>
    eqn.B1 = rand(500, 2);

    rand('seed', 6.0); %#ok<RAND>
    eqn.C2 = rand(3, 500);
    rand('seed', 7.0); %#ok<RAND>
    eqn.C1 = rand(3, 500);
else
    rng(1.0);
    A = rand(500) - 250 * eye(500);
    rng(2.0);
    M = rand(500);
    M = M' * M;
    rng(3.0);
    J = rand(100, 500);

    rng(4.0);
    eqn.B2 = rand(500, 2);
    rng(5.0);
    eqn.B1 = rand(500, 2);

    rng(6.0);
    eqn.C2 = rand(3, 500);
    rng(7.0);
    eqn.C1 = rand(3, 500);
end

eqn.A_ = sparse([A, J'; J, zeros(100)]);
eqn.E_ = sparse(blkdiag(M, zeros(100)));
st     = 500;
eqn.st = 500;

gam    = 5;
eqn.B1 = 1/gam * eqn.B1;

eqn.type  = 'T';
eqn.haveE = 1;

%% Set operator.
oper = operatormanager('dae_2');

%% Construction of options struct.
% ADI settings.
opts.adi.maxiter          = 200;
opts.adi.res_tol           = 1.0e-14;
opts.adi.rel_diff_tol            = 0;
opts.adi.info             = 1;
opts.adi.compute_sol_fac  = 1;
opts.adi.accumulateK      = 0;
opts.adi.accumulateDeltaK = 0;

% Shift options.
opts.shifts.num_desired     = 5;
opts.shifts.method = 'projection';

% % NM settings.
opts.nm.maxiter       = 50;
opts.nm.res_tol        = 1.0e-12;
opts.nm.rel_diff_tol         = 1.0e-12;
opts.nm.info          = 1;
opts.nm.linesearch    = 0;
opts.nm.accumulateRes = 1;

% NM projection settings.
opts.nm.projection      = [];
opts.nm.projection.freq = 0;
opts.nm.res.maxiter     = 10;
opts.nm.res.tol         = 1.0e-06;
opts.nm.res.orth        = 1;

% RI settings.
opts.ri.riccati_solver = 'newton';
opts.ri.maxiter        = 10;
opts.ri.res_tol         = 1.0e-09;
opts.ri.rel_diff_tol          = 1.0e-16;
opts.ri.compres_tol    = 1.0e-16;
opts.ri.info           = 1;

% global options
opts.norm           = 2;

% %% Call Riccati iteration with Newton solver.
t_mess_lrri = tic;
[outnm, eqn, opts, oper] = mess_lrri(eqn, opts, oper);
t_elapsed1 = toc(t_mess_lrri);
fprintf(1,'mess_lrri took %6.2f seconds \n' , t_elapsed1);
%% Setup RADI structure.
opts.radi.maxiter = opts.adi.maxiter;
opts.radi.res_tol  = opts.nm.res_tol;
opts.radi.rel_diff_tol   = 1.0e-16;
opts.radi.info    = 1;

opts.ri.riccati_solver = 'radi';

%% Call Riccati iteration with RADI solver.
t_mess_lrri = tic;
[out, eqn, opts, ~] = mess_lrri(eqn, opts, oper);
t_elapsed2 = toc(t_mess_lrri);
fprintf(1,'mess_lrri took %6.2f seconds \n',t_elapsed2);
%% Test of the solution.
% Partitioning of the system.
A = eqn.A_(1:st,1:st);
J = eqn.A_(1:st,st+1:end);
G = eqn.A_(st+1:end,1:st);
E = eqn.E_(1:st,1:st);
B1 = eqn.B1;
B2 = eqn.B2;
C1 = eqn.C1;

% Compute projection matrices (not recommended for large-scale case).
Pi_l = eye(st) - J*((G*(E\J))\(G/E));
Pi_r = eye(st) - (E\J)*((G*(E\J))\G);

% Explicit projection.
A_p  = Pi_l * A * Pi_r;
M_p  = Pi_l * E * Pi_r;
C1_p = C1 * Pi_r;
B1_p = Pi_l * B1;
B2_p = Pi_l * B2;

% Compute the actual errors.
abserrnm = norm(A_p' * (outnm.Z * outnm.Z') * M_p ...
    + M_p' * (outnm.Z * outnm.Z') * A_p ...
    + M_p' * (outnm.Z * outnm.Z') * (B1_p * B1_p' ...
    - B2_p * B2_p') * (outnm.Z * outnm.Z') * M_p + C1_p' * C1_p, 2);
relerrnm = abserrnm / norm(C1_p * C1_p', 2);
fprintf(1, '\nNewton -> set tolerance vs. real residual: %e | %e\n', ...
    opts.ri.res_tol, relerrnm);

abserrradi = norm(A_p' * (out.Z * out.Z') * M_p ...
    + M_p' * (out.Z * out.Z') * A_p ...
    + M_p' * (out.Z * out.Z') * (B1_p * B1_p' ...
    - B2_p * B2_p') * (out.Z * out.Z') * M_p + C1_p' * C1_p, 2);
relerrradi = abserrradi / norm(C1_p * C1_p', 2);
fprintf(1, 'RADI   -> set tolerance vs. real residual: %e | %e\n', ...
    opts.ri.res_tol, relerrradi);

if istest
    assert(relerrnm < opts.ri.res_tol, ...
        'MESS:TEST:accuracy','unexpectedly inaccurate result');
    assert(relerrradi < opts.ri.res_tol, ...
        'MESS:TEST:accuracy','unexpectedly inaccurate result');
end
