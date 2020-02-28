function DEMO_RI_GE_T_N(istest)
% Computes the solution of the Hinf Riccati equation for a random generated
% generalized system. The computations are done with the RADI method for
% the control and filter Hinf Riccati equations. Afterwards, the real
% residual norms are shown and compared to the set tolerance.
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
    eqn.E_ = rand(500);
    eqn.E_ = eqn.E_' * eqn.E_;
    
    rand('seed', 3.0); %#ok<RAND>
    eqn.B2 = rand(500, 2);
    rand('seed', 4.0); %#ok<RAND>
    B1 = rand(500, 2);
    
    rand('seed', 5.0); %#ok<RAND>
    eqn.C2 = rand(3, 500);
    rand('seed', 6.0); %#ok<RAND>
    C1 = rand(3, 500);
else
    rng(1.0);
    eqn.A_ = rand(500) - 250 * eye(500);
    rng(2.0);
    eqn.E_ = rand(500);
    eqn.E_ = eqn.E_' * eqn.E_;
    
    rng(3.0);
    eqn.B2 = rand(500, 2);
    rng(4.0);
    B1 = rand(500, 2);
    
    rng(5.0);
    eqn.C2 = rand(3, 500);
    rng(6.0);
    C1 = rand(3, 500);
end

eqn.haveE = 1;

gam = 5; % Scaling term for disturbances.

%% Set operator.
oper = operatormanager('default');

%% Construction of options struct.
% RADI settings.
opts.radi.maxiter      = 100;
opts.radi.res_tol      = 1.0e-12;
opts.radi.rel_diff_tol = 1.0e-16;
opts.radi.info         = 1;
opts.radi.trunc_tol    = eps;

% Shift options.
opts.shifts.num_desired     = 5;
opts.shifts.method = 'projection';

% RI settings.
opts.ri.riccati_solver = 'radi';
opts.ri.maxiter        = 10;
opts.ri.res_tol        = 1.0e-09;
opts.ri.rel_diff_tol   = 1.0e-16;
opts.ri.compres_tol    = 1.0e-16;
opts.ri.info           = 1;
opts.ri.trunc_tol      = eps;

% global options
opts.norm           = 2;

%% Solve the control equation.
tic;
eqn.type = 'T';
eqn.B1   = 1/gam * B1;
eqn.C1   = C1;
[outControl, eqn, opts, oper] = mess_lrri(eqn, opts, oper);
toc;

%% Solve the filter equation.
tic;
eqn.type = 'N';
eqn.B1   = B1;
eqn.C1   = 1/gam * C1;
[outFilter, eqn, opts, ~] = mess_lrri(eqn, opts, oper);
toc;

%% Compute real residuals.
absControl = norm(eqn.A_' * (outControl.Z * outControl.Z') * eqn.E_ ...
    + eqn.E_' * (outControl.Z * outControl.Z') * eqn.A_ ...
    + eqn.E_' * (outControl.Z * outControl.Z') * (1/gam^2 * (B1 * B1') ...
    - eqn.B2 * eqn.B2') * (outControl.Z * outControl.Z') * eqn.E_ ...
    + C1' * C1, 2);
relControl = absControl / norm(C1 * C1', 2);
fprintf(1, '\nControl -> set tolerance vs. real residual: %e | %e\n', ...
    opts.ri.res_tol, relControl);

absFilter = norm(eqn.A_ * (outFilter.Z * outFilter.Z') * eqn.E_' ...
    + eqn.E_ * (outFilter.Z * outFilter.Z') * eqn.A_' ...
    + eqn.E_ * (outFilter.Z * outFilter.Z') * (1/gam^2 * (C1' * C1) ...
    - eqn.C2' * eqn.C2) * (outFilter.Z * outFilter.Z') * eqn.E_' ...
    + B1 * B1', 2);
relFilter = absFilter / norm(B1' * B1, 2);
fprintf(1, 'Filter  -> set tolerance vs. real residual: %e | %e\n', ...
    opts.ri.res_tol, relFilter);

if istest
    assert(relControl < opts.ri.res_tol, ...
        'MESS:TEST:accuracy','unexpectedly inaccurate result');
    assert(relFilter < opts.ri.res_tol, ...
        'MESS:TEST:accuracy','unexpectedly inaccurate result');
end
