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
%               2009-2020
%
%%
if nargin<1, istest=0; end

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
oper = operatormanager('default');

%% Construction of options struct.
% ADI settings.
opts.adi.maxiter          = 200;
opts.adi.res_tol           = 1.0e-12;
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
opts.nm.res_tol        = 1.0e-10;
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

% RADI settings.
opts.radi.maxiter = opts.adi.maxiter;
opts.radi.res_tol  = opts.nm.res_tol;
opts.radi.rel_diff_tol   = 1.0e-16;
opts.radi.info    = 1;

% RI settings.
opts.ri.riccati_solver = 'radi';
opts.ri.lqg_solver     = 'newton';
opts.ri.maxiter        = 10;
opts.ri.res_tol         = 1.0e-10;
opts.ri.rel_diff_tol          = 1.0e-16;
opts.ri.compres_tol    = 1.0e-16;
opts.ri.info           = 1;

% global options
opts.norm           = 2;

% %% Call Riccati iteration with Newton solver.
tic;
eqn.type = 'T';
eqn.B1   = 1/gam * B1;
eqn.C1   = C1;
[out, eqn, opts, ~] = mess_lrri(eqn, opts, oper);
toc;

%% Compute real residuals.
abserr = norm(eqn.A_' * (out.Z * out.Z') + (out.Z * out.Z') * eqn.A_ ...
    + (out.Z * out.Z') * (1/gam^2 * (B1 * B1') ...
    - eqn.B2 * eqn.B2') * (out.Z * out.Z') + C1' * C1, 2);
relerr = abserr / norm(C1 * C1', 2);
fprintf(1, '\nset tolerance vs. real residual: %e | %e\n', ...
    opts.ri.res_tol, relerr);

if istest
    assert(relerr < opts.ri.res_tol, ...
        'MESS:TEST:accuracy','unexpectedly inaccurate result');
end
