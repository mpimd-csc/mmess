function IRKA_mor_Stokes(istest)
% Computes a standard ROM by implicitly running IRKA for the equivalent
% projected system on the hidden manifold.
%
% Inputs:
% istest        flag to determine whether this demo runs as a CI test or
%               interactive demo
%               (optional, defaults to 0, i.e. interactive demo)
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% IRKA tolerance and maximum iteration number
opts.irka.maxiter = 150;
opts.irka.r = 30;
opts.irka.flipeig = false;
opts.irka.h2_tol = 1e-6;

opts.irka.init = 'logspace';

if nargin < 1
    istest = false;
end

if istest
    opts.irka.info = 1;
else
    opts.irka.info = 2;
end

[oper, opts] = operatormanager(opts, 'dae_2');

%% Problem data
nin = 5;
nout = 5;
nx = 10;
ny = 10;
[eqn.E_, eqn.A_, eqn.Borig, eqn.Corig] = ...
    stokes_ind2(nin, nout, nx, ny, opts);
n = size(eqn.E_, 1);
eqn.haveE = true;
n_ode = trace(eqn.E_); % Stokes is FDM discretized, so so this is
% the dimension of the velocity space
eqn.manifold_dim = n_ode;
eqn.B = eqn.Borig(1:n_ode, :);
eqn.C = eqn.Corig(:, 1:n_ode);
%% Compute reduced system matrices
t_mess_tangential_irka = tic;
[ROM.E, ROM.A, ROM.B, ROM.C, ~, outinfo] = ...
    mess_tangential_irka(eqn, opts, oper);
t_elapsed1 = toc(t_mess_tangential_irka);
mess_fprintf(opts, 'mess_tangential_irka took %6.2f seconds \n\n', t_elapsed1);

%%
t_eval_ROM = tic;

%% Evaluate the ROM quality
% while the Gramians are computed exploiting the DAE structure, due to the
% construction of the function handles we can not do so for the transfer
% function. Therefore we need to extend the matrices B and C and call the
% 'default' usfs for unstructured computation:
eqn.B = eqn.Borig;
eqn.C = eqn.Corig;
[oper, opts] = operatormanager(opts, 'default');

if istest
    opts.tf_plot.info = 0;
else
    opts.tf_plot.info = 2;
end

opts.tf_plot.fmin = -3;
opts.tf_plot.fmax = 4;

opts.tf_plot.type = 'sigma';

% We want to have the frequencies depicted in Hertz and values in decibels
opts.tf_plot.Hz = true;
opts.tf_plot.db = true;

out = mess_tf_plot(eqn, opts, oper, ROM);

t_elapsed2 = toc(t_eval_ROM);
mess_fprintf(opts, ...
             'evaluation of ROM matrices took %6.2f seconds \n', t_elapsed2);

%%
maerr = max(abs(out.err));
mrerr = max(abs(out.relerr));
if istest && (maerr > 1e-6 || mrerr > 1e-4)
    mess_err(opts, 'TEST:accuracy', ['unexpectedly inaccurate result.\n' ...
                                     'max. abs err: %e (allowed 1e-6)\n' ...
                                     'max rel err: %e (allowed 1e-4)'], ...
             maerr, mrerr);
end

%%
problem = 'Stokes';
fprintf(['\nComputing open loop step response of original and ', ...
         'reduced-order systems and time domain MOR errors\n']);
open_step(eqn, ROM.A, ROM.B, ROM.C, problem, istest);

%%
fprintf('\nComputing ROM based feedback\n');
if exist('care', 'file')
    [~, ~, Kr] = care(ROM.A, ROM.B, ROM.C' * ROM.C, eye(size(ROM.B, 2)));
else
    Y = care_nwt_fac([], ROM.A, ROM.B, ROM.C, 1e-12, 50);
    Kr = (Y * ROM.B)' * Y;
end
K = [Kr * outinfo.TL' * eqn.E_(1:n_ode, 1:n_ode), zeros(size(Kr, 1), n - n_ode)];

%%
fprintf(['\nComputing closed loop step response of original and ', ...
         'reduced-order systems and time domain MOR errors\n']);
closed_step(eqn, ROM.A, ROM.B, ROM.C, problem, K, Kr, istest);
