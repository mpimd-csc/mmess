function LQR_TripleChain(n1, usfs, shifts, istest)
% LQR_TripleChain computes the optimal feedback for the LQR problem (e.g.
% [1]) with respect to the Triple-Chain-Oscillator from [2]. The
% computations use either first or second companion form linearization. In
% order to reduce the complexity, the 'so_1' and 'so_2' function handles
% cast all operations in the linearized 2n x 2n model back to operations
% with respect to the original n x n matrices [3,4]. The default function
% handles explicitly form the 2n x 2n matrices and treat them as a first
% order system.
%
% Usage:      LQR_TripleChain(n1, oper, shifts, istest)
%
% Inputs:
%
% n1          length of a single chain in the model
%             (optional; defaults to 1000)
%
% usfs        the set of user supplied functions to use.
%             possible values: 'so_1', 'so_2', 'default'
%             (optional; defaults to 'so_1'
%
% shifts      the desired ADI shift selection strategy
%             possible values: 'heur', 'projection'
%             (optional; defaults to 'projection'
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% References:
%
% [1] A. Locatelli, Optimal Control: An Introduction, Birkhäuser, Basel,
%     Switzerland, 2001.
%
% [2] N. Truhar, K. Veselić, An efficient method for estimating the optimal
%     dampers’ viscosity for linear vibrating systems using Lyapunov
%     equation, SIAM J. Matrix Anal. Appl. 31 (1) (2009) 18–39.
%     https://doi.org/10.1137/070683052
%
% [3] P. Benner, J. Saak, Efficient Balancing based MOR for Second Order
%     Systems Arising in Control of Machine Tools, in: I. Troch, F.
%     Breitenecker (Eds.), Proceedings of the MathMod 2009, no. 35 in
%     ARGESIM-Reports, Vienna Univ. of Technology, ARGE Simulation News,
%     Vienna, Austria, 2009, pp. 1232–1243, iSBN/ISSN:978-3-901608-35-3
%
% [4] P. Benner, P. Kürschner, J. Saak, An improved numerical method for
%     balanced truncation for symmetric second order systems,
%     Mathematical and Computer Modelling of Dynamical Systems 19 (6)
%     (2013) 593–615.
%     https://doi.org/10.1080/13873954.2013.794363

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
narginchk(0, 4);
if nargin < 1
    n1 = 1000;
end
if nargin < 2
    usfs = 'so_1';
end
if nargin < 3
    shifts = 'projection';
end
if nargin < 4
    istest = false;
end
%%
% set operation
opts = struct();
[oper, opts] = operatormanager(opts, usfs);
% Problem data

alpha = 2;
Beta  = 5;
v     = 5;

switch usfs
    case 'so_1'
        [eqn.M_, eqn.E_, eqn.K_] = ...
            triplechain_MSD(n1, alpha, Beta, v);
        s  = size(eqn.K_, 1);
        eqn.B = [zeros(s, 1); ones(size(eqn.K_, 1), 1)];
        eqn.C = [ones(1, size(eqn.K_, 1)) zeros(1, s)];
    case 'so_2'
        [eqn.M_, eqn.E_, eqn.K_] = ...
            triplechain_MSD(n1, alpha, Beta, v);
        s  = size(eqn.K_, 1);
        eqn.B = [ones(size(eqn.K_, 1), 1); zeros(s, 1)];
        eqn.C = eqn.B';
    case 'default'
        [M_, E_, K_] = triplechain_MSD(n1, alpha, Beta, v);
        s = size(K_, 1);
        eqn.A_ = [sparse(s, s), -K_; -K_, -E_];
        eqn.E_ = [-K_, sparse(s, s); sparse(s, s), M_];
        eqn.B = [zeros(s, 1); ones(size(K_, 1), 1)];
        eqn.C = [ones(1, size(K_, 1)) zeros(1, s)];
        clear M_ E_ K_ s;
end

eqn.haveE = true;

%%
% ADI tolerances and maximum iteration number
opts.adi.maxiter          = 200;
opts.adi.res_tol          = 1e-10;
opts.adi.rel_diff_tol     = 0;
opts.adi.info             = 0;
opts.adi.accumulateK      = true;
opts.adi.accumulateDeltaK = false;
opts.adi.compute_sol_fac  = true;

eqn.type = 'T';
%%
% Heuristic shift parameters via basic Arnoldi
n = oper.size(eqn, opts);
switch shifts
    case 'heur'
        opts.shifts.num_desired = 25;
        opts.shifts.num_Ritz    = 50;
        opts.shifts.num_hRitz   = 25;

        opts.shifts.info        = 0;
        opts.shifts.method      = 'heur';
        opts.shifts.b0          = ones(n, 1);
    case 'projection'
        opts.shifts.num_desired = 6;
        opts.shifts.method      = 'projection';
        n                       = oper.size(eqn, opts);
        opts.shifts.b0          = ones(n, 1);
end
opts.shifts.truncate = 1e6; % remove all shifts larger than 1e6 or smaller
% than 1e-6 in absolute value in order to avoid
% loosing information about M or K in the
% shifted coefficients (p^2*M-pD+K)
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter                 = 25;
opts.nm.res_tol                 = 1e-10;
opts.nm.rel_diff_tol            = 1e-16;
opts.nm.info                    = 1;
opts.nm.accumulateRes           = true;
opts.nm.linesearch              = true;
opts.nm.inexact                 = 'quadratic';
opts.nm.tau                     = 0.1;
opts.norm                       = 'fro';

%%
t_mess_lrnm = tic;
outnm = mess_lrnm(eqn, opts, oper);
t_elapsed1 = toc(t_mess_lrnm);
mess_fprintf(opts, 'mess_lrnm took %6.2f seconds \n', t_elapsed1);
if istest
    if min(outnm.res) >= opts.nm.res_tol
        mess_err(opts, 'TEST:accuracy', 'unexpectedly inaccurate result');
    end
else
    figure(1);
    semilogy(outnm.res, 'LineWidth', 3);
    title('0= C^TC + A^T X E + E^T X A -E^T X BB^T X M');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
[mZ, nZ] = size(outnm.Z);
mess_fprintf(opts, 'size outnm.Z: %d x %d\n\n', mZ, nZ);

%% Lets try the RADI method and compare
% RADI-MESS settings
opts.shifts.history = opts.shifts.num_desired * size(eqn.C, 1);

opts.shifts.method = 'projection';

% .. Suggest false (smart update is faster; convergence is the same).
opts.shifts.naive_update_mode  = false;
opts.radi.compute_sol_fac      = true;
opts.radi.get_ZZt              = true;
opts.radi.maxiter              = opts.adi.maxiter;
opts.norm                      = 2;
opts.radi.res_tol              = opts.nm.res_tol;
opts.radi.rel_diff_tol         = 0;
opts.radi.info                 = 1;

t_mess_lrradi = tic;
outradi = mess_lrradi(eqn, opts, oper);
t_elapsed2 = toc(t_mess_lrradi);
mess_fprintf(opts, 'mess_lrradi took %6.2f seconds \n', t_elapsed2);
if istest
    if min(outradi.res) >= opts.radi.res_tol
        mess_err(opts, 'TEST:accuracy', 'unexpectedly inaccurate result');
    end
else
    figure(2);
    semilogy(outradi.res, 'LineWidth', 3);
    title('0= C^T C + A^T X E + E^T X A -E^T X BB^T X E');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
end
[mZ, nZ] = size(outradi.Z);
mess_fprintf(opts, 'size outradi.Z: %d x %d\n\n', mZ, nZ);

%% compare
if istest
    nrm = norm(outnm.K - outradi.K, 'fro');
    nrmNM = norm(outnm.K, 'fro');
    if nrm / nrmNM >= 1e-9
        mess_err(opts, 'TEST:accuracy', ...
                 'unexpectedly inaccurate result: ', ...
                 '||K_NM - K_RADI||_F / ||K_NM||_F=%g', nrm / nrmNM);
    end
else
    figure(3);
    ls_nm   = [outnm.adi.niter];
    ls_radi = 1:outradi.niter;

    semilogy(cumsum(ls_nm), outnm.res, 'k--', ...
             ls_radi, outradi.res, 'b-', ...
             'LineWidth', 3);
    title('0= C^T C + A^T X E + E^T X A - E^T X BB^T X E');
    xlabel('number of solves with A+p*E');
    ylabel('normalized residual norm');
    legend('LR-NM', 'RADI');
end
