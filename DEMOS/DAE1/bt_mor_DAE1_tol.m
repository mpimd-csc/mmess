function bt_mor_DAE1_tol(k, istest)
% Computes a reduced order model via Balanced Truncation for the proper
% index-1 System BIPS98_606 from https://sites.google.com/site/rommes/software
% following the method suggested in [1]
%
% Input:
% k         select model (allowed values 1,..5, 7, .., 13, default 7)
%           note that BIPS model 6 is not stable.
% istest    decides whether the function runs as an interactive demo or a
%           continuous integration test. (optional; defaults to 0, i.e.
%           interactive demo)
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
opts = struct;

%%
if nargin < 1
    k = 7;
end
if nargin < 2
    istest = false;
end
%%
% set operation manager for the Gramian computations
[oper, opts] = operatormanager(opts, 'dae_1');

%% Read problem data
if k == 6
    mess_warn(opts, 'illegal_input', ...
              ['The Juba5723 model is not stable and thus not supported. ', ...
               'Using BIPS bips98_606 instead.']);
    k = 7;
end
eqn = mess_get_BIPS(k);

%% Turn off  close to singular warnings
% (this model is really badly conditioned)
orig_warnstate = warning('OFF', 'MATLAB:nearlySingularMatrix');

%%
% ADI tolerances and maximum iteration number
opts.adi.maxiter       = 300;
opts.adi.res_tol       = 1e-10;
opts.adi.rel_diff_tol  = 1e-16;
opts.adi.info          = 1;
opts.norm = 'fro';

%%
opts.shifts.method = 'projection';
opts.shifts.num_desired = 20;

%%
eqn.type = 'N';
t_mess_lradi = tic;
outB = mess_lradi(eqn, opts, oper);
t_elapsed1 = toc(t_mess_lradi);
mess_fprintf(opts, 'mess_lradi took %6.2f seconds \n', t_elapsed1);
if istest
    if min(outB.res) >= opts.adi.res_tol
        mess_err(opts, 'TEST:accuracy', 'unexpectedly inaccurate result');
    end
else
    figure;
    semilogy(outB.res, 'LineWidth', 3);
    title('0 = BB^T + A X E^T + E X A^T');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
end

[mZ, nZ] = size(outB.Z);
mess_fprintf(opts, 'size outB.Z: %d x %d\n', mZ, nZ);

%%
eqn.type     = 'T';
t_mess_lradi = tic;
outC         = mess_lradi(eqn, opts, oper);
t_elapsed2   = toc(t_mess_lradi);

mess_fprintf(opts, 'mess_lradi took %6.2f seconds \n', t_elapsed2);

if istest
    if min(outC.res) >= opts.adi.res_tol
        mess_err(opts, 'TEST:accuracy', ...
                 'unexpectedly inaccurate result');
    end
else
    figure;
    semilogy(outC.res, 'LineWidth', 3);
    title('0 = C^T C + A^T X E + E^T X A');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
[mZ, nZ] = size(outC.Z);
mess_fprintf(opts, 'size outC.Z: %d x %d\n', mZ, nZ);

%% Compute reduced system matrices
% Perform Square Root Method  (SRM)

% BT tolerance and maximum order for the ROM
opts.srm.tol     = 1e-3;
opts.srm.max_ord = 250;

% SRM verbosity
if istest
    opts.srm.info = 1;
else
    opts.srm.info = 2;
end

% The actual SRM
[TL, TR, hsv] = mess_square_root_method(eqn, opts, oper, outB.Z, outC.Z);

% compute ROM matrices
n = size(eqn.A_, 1);
one = 1:eqn.manifold_dim;
two = (eqn.manifold_dim + 1):n;

B1 = TL' * (eqn.A_(one, one) * TR);
B2 = TL' * (eqn.A_(one, two));
A1 = eqn.A_(two, one) * TR;

ROM.A = B1 - B2 * (eqn.A_(two, two) \ A1);
ROM.B = TL' * eqn.B(one, :) - B2 * (eqn.A_(two, two) \ eqn.B(two, :));
ROM.C = eqn.C(:, one) * TR - eqn.C(:, two) * (eqn.A_(two, two) \ A1);
ROM.D = -eqn.C(:, two) * (eqn.A_(two, two) \ eqn.B(two, :));
ROM.E = eye(size(ROM.A));

%% Evaluate the ROM quality
% while the Gramians are computed on the hidden manifold, we need to do the
% frequency domain computations without (implicitly) using the Schur
% complement (due to the construction of the function handles)
[oper, opts] = operatormanager(opts, 'default');

if istest
    opts.tf_plot.info = 0;
else
    opts.tf_plot.info = 2;
end

opts.tf_plot.fmin = -3;
opts.tf_plot.fmax = 4;

opts.tf_plot.type = 'sigma';

out = mess_tf_plot(eqn, opts, oper, ROM);
err = out.err;

if istest
    if max(err) > 5e-3
        mess_err(opts, 'TEST:accuracy', 'unexpectedly inaccurate result');
    end
else
    figure;
    semilogy(hsv, 'LineWidth', 3);
    title('Computed Hankel singular values');
    xlabel('index');
    ylabel('magnitude');
end
%% reset warning state
warning(orig_warnstate);
