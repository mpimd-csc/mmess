function LQR_DAE1(istest)
% Computes a Riccati feedback control for the proper
% index-1 System BIPS98_606 from https://sites.google.com/site/rommes/software
% following the ideas introduced in [1] for Lyapunov equations using he
% Newton-ADI iteration.
%
% Input:
% istest  decides whether the function runs as an interactive demo or a
%         continuous integration test. (optional; defaults to 0, i.e.
%         interactive demo)
%
% References:
% [1] F. Freitas, J. Rommes, N. Martins, Gramian-based reduction method
%    applied to large sparse power system descriptor models, IEEE Transactions
%    on Power Systems 23 (3) (2008) 1258–1270. doi:10.1109/TPWRS.2008.926693

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

%%
% set operation
opts = struct();
[oper, opts] = operatormanager(opts, 'dae_1');

%% Problem data
eqn = mess_get_BIPS(7);

%% Turn off close to singular warnings
%  (this model is really badly conditioned)
orig_warnstate = warning('OFF', 'MATLAB:nearlySingularMatrix');

%%
opts.norm = 'fro';

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 300;
opts.adi.res_tol = 1e-12;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;
opts.adi.projection.freq = 0;

eqn.type = 'N';
%%
opts.shifts.num_desired = 25;
opts.shifts.num_Ritz = 50;
opts.shifts.num_hRitz = 25;
opts.shifts.method = 'projection';
opts.shifts.num_desired = 9;

%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 30;
opts.nm.res_tol = 1e-10;
opts.nm.rel_diff_tol = 1e-16;
opts.nm.info = 1;
opts.nm.linesearch = true;
opts.nm.accumulateRes = true;
%%
t_mess_lrnm = tic;
outnm = mess_lrnm(eqn, opts, oper);
t_elapsed1 = toc(t_mess_lrnm);
mess_fprintf(opts, 'mess_lrnm took %6.2f seconds \n', t_elapsed1);
if istest
    if min(outnm.res) >= opts.nm.res_tol
        mess_err(opts, 'TEST:accuracy', ...
                 'unexpectedly inaccurate result in LRNM');
    end
else
    figure();
    semilogy(outnm.res, 'LineWidth', 3);
    title('0 = C^T C + A^T X E + E^T X A -E^T X BB^T X E');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
[mZ, nZ] = size(outnm.Z);
mess_fprintf(opts, 'size outnm.Z: %d x %d\n', mZ, nZ);

%% Lets try RADI
opts.norm = 2;

% RADI-MESS settings
opts.shifts.history = opts.shifts.num_desired * size(eqn.C, 1);
opts.shifts.method  = 'gen-ham-opti';

opts.shifts.naive_update_mode = false;
% .. Suggest false (smart update is faster; convergence is the same).

opts.radi.compute_sol_fac = true;
opts.radi.get_ZZt         = true;
opts.radi.compute_res     = false;
opts.radi.maxiter         = 500;
opts.radi.res_tol         = opts.nm.res_tol;
opts.radi.rel_diff_tol    = 0;
opts.radi.info            = 1;

t_mess_lrradi = tic;
outradi = mess_lrradi(eqn, opts, oper);
t_elapsed2 = toc(t_mess_lrradi);
mess_fprintf(opts, 'mess_lrradi took %6.2f seconds \n', t_elapsed2);
if istest
    if min(outradi.res) >= opts.radi.res_tol
        mess_err(opts, 'TEST:accuracy', ...
                 'unexpectedly inaccurate result in RADI');
    end
else
    figure();
    semilogy(outradi.res, 'LineWidth', 3);
    title('0 = C^T C + A^T X E + E^T X A - E^T X BB^T X E');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
end
[mZ, nZ] = size(outradi.Z);
mess_fprintf(opts, 'size outradi.Z: %d x %d\n', mZ, nZ);

%% compare
if not(istest)
    figure();
    ls_nm = [outnm.adi.niter];
    ls_radi = 1:outradi.niter;

    semilogy(cumsum(ls_nm), outnm.res, 'k--', ...
             ls_radi, outradi.res, 'b-', ...
             'LineWidth', 3);
    title('0 = C^T C + A^T X E + E^T X A -E^T X BB^T X E');
    xlabel('number of solves with A + p * E');
    ylabel('normalized residual norm');
    legend('LR-NM', 'RADI');
end

%% reset warning state
warning(orig_warnstate');
