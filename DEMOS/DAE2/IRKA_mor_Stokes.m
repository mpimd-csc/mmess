function IRKA_mor_Stokes(istest)
% Computes a standard ROM by implicitly running IRKA for the equivalent
% projected system on the hidden manifold.
%
% Inputs:
% istest        flag to determine whether this demo runs as a CI test or
%               interactive demo
%               (optional, defaults to 0, i.e. interactive demo)
%


% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others
%               2009-2019
%

% IRKA tolerance and maximum iteration number
opts.irka.maxiter = 150;
opts.irka.r = 30;
opts.irka.flipeig = 0;
opts.irka.h2_tol = 1e-6;

opts.irka.init = 'logspace';

if nargin < 1, istest = 0; end
if istest
    opts.irka.info = 1;
else
    opts.irka.info = 1;
end

oper = operatormanager('dae_2');

%% Problem data
nin = 5;
nout = 5;
nx = 10;
ny = 10;
[eqn.E_, eqn.A_, eqn.Borig, eqn.Corig] = ...
    stokes_ind2(nin, nout, nx, ny);
n = size(eqn.E_, 1);
eqn.haveE = 1;
st = full(sum(diag(eqn.E_)));
eqn.st = st;
eqn.B = eqn.Borig(1:st, :);
eqn.C = eqn.Corig(:, 1:st);
%% Compute reduced system matrices
tic;
[ROM.E, ROM.A, ROM.B, ROM.C, S, b, c, V, W] = ...
    mess_tangential_irka(eqn, opts, oper);
toc;

%%
tic;

%% Evaluate the ROM quality
% while the Gramians are computed exploiting the DAE structure, due to the
% construction of the function handles we can not do so for the transfer
% function. Therfore we need to extend the matrices B and C and call the
% 'default' usfs for unstructured computation:
eqn.B = eqn.Borig;
eqn.C = eqn.Corig;
oper = operatormanager('default');

if istest
    opts.sigma.info = 0;
else
    opts.sigma.info = 2;
end

opts.sigma.fmin = -3;
opts.sigma.fmax = 4;

out = mess_sigma_plot(eqn, opts, oper, ROM);

toc;
%%
maerr = max(abs(out.err));
mrerr = max(abs(out.relerr));
if istest && (maerr>1e-6 || mrerr>1e-4)
    error('MESS:TEST:accuracy',['unexpectedly inaccurate result.\n' ...
                        'max. abs err: %e (allowed 1e-6)\n' ...
                        'max rel err: %e (allowed 1e-4)'], maerr, mrerr);
end

%%
problem = 'Stokes';
fprintf(['\nComputing open loop step response of original and ', ...
    'reduced-order systems and time domain MOR errors\n']);
open_step(eqn, ROM.A, ROM.B, ROM.C, problem, istest);

%%
fprintf('\nComputing ROM based feedback\n');
if exist('care', 'file')
    [~, ~, Kr] = care(ROM.A, ROM.B, ROM.C'*ROM.C, eye(size(ROM.B, 2)));
else
    Y = care_nwt_fac([], ROM.A, ROM.B, ROM.C, 1e-12, 50);
    Kr = (Y * ROM.B)' * Y;
end
K = [Kr * W' * eqn.E_(1:st, 1:st), zeros(size(Kr, 1), n-st)];

%%
fprintf(['\nComputing closed loop step response of original and ', ...
    'reduced-order systems and time domain MOR errors\n']);
closed_step(eqn, ROM.A, ROM.B, ROM.C, problem, K, Kr, istest);
