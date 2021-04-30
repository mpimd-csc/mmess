function LQG_FDM_unstable_radi(n0, n_unstable, istest)
% Computes stabilizing and detecting solutions of the regulator and filter
% Riccati equations via the RADI, respectively. The solution of the
% regulator equation is approximated via a ZZ' factorization, and for the
% filter equation an LDL' approximation is used.
% Note: The LDL' case is only used for demonstrational reasons but with no
% underlying necessity here.
%
% Inputs:
%
% n0          n0^2 gives the dimension of the original model, i.e. n0 is
%             the number of degrees of freedom, i.e. grid points, per
%             spatial direction
%             (optional; defaults to 30)
%
% n_unstable  number of unstable eigenvalues by construction
%             (optional; defaults to 5)
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

narginchk(0,3);
if nargin<1, n0 = 30; end
if nargin<2, n_unstable = 5; end
if nargin<3, istest = 0; end


%% construct unstable matrix
% Generate the stable part.
A0         = fdm_2d_matrix(n0, '10*x', '1000*y', '0');
n_stable   = size(A0, 1);
n          = n_stable + n_unstable;
B          = [eye(n_stable, n_unstable); diag(1:n_unstable)];

C                                 = [transpose(B); zeros(2, n)];
C(n_unstable + 1, n_unstable + 1) = 1;
C(n_unstable + 2, n_unstable + 2) = 1;

% The instability is added by fixing the solution Y of the partial
% stabilization problem via Bernoulli and solving the resulting Lyapunov
% equation for the unstable matrix Au:
%
%   Au'*Y + Y*Au - Y*Bt*Bt'*Y = 0.
%
% The resulting stabilizing feedback is defined such that the unstable
% eigenvalues of the system are mirrored into the left open half-plane.

Y   = eye(n_unstable);
Bt  = B(n_stable+1:n, : ); % Part of B acting on unstable part.
Vl  = [zeros(n_stable, n_unstable); eye(n_unstable)];
if exist('lyap','file') % Solve Lyapunov equation.
    Au = lyap(Y, -Bt * Bt');
else
    Au = lyap2solve(Y, (-Bt * Bt')');
end

ZB0 = Vl;
YB0 = Y \ eye(n_unstable);
WB0 = C';

ZC0 = ZB0;
YC0 = YB0;
WC0 = B;

% The full A is constructed via additive decomposition (blockdiagonal).
eqn.A_     = blkdiag(A0,Au);
eqn.B      = B;
eqn.C      = C;
eqn.haveE  = 0;

% set operation
oper = operatormanager('default');


%% global options
opts.norm = 'fro';


%% RADI parameters
opts.shifts.naive_update_mode = 0;
opts.radi.compute_sol_fac     = 1;
opts.radi.compute_res         = 0;
opts.radi.get_ZZt             = 1;
opts.radi.maxiter             = 200;
opts.radi.res_tol             = 1.0e-10;
opts.radi.rel_diff_tol        = 1.0e-16;
opts.radi.info                = 1;

% Only for mess_res2_norms (not needed in RADI):
opts.radi.res = struct('maxiter', 10, 'tol', 1.0e-06, 'orth', 0);


%% shift parameters via projection
opts.shifts.num_desired = 5;
opts.shifts.method      = 'projection';


%% for shift banned eigenvalues
opts.shifts.banned     = -eig(Au);
opts.shifts.banned_tol = 1e-6;


%% Solve regulator Riccati equation.
fprintf('Solve regulator Riccati equation\n');
fprintf('--------------------------------\n');
eqn.type     = 'T';
opts.LDL_T   = 0;
opts.radi.Z0 = ZB0;
opts.radi.Y0 = YB0;
opts.radi.W0 = WB0;

tic;
[outReg, eqn, opts, oper] = mess_lrradi(eqn, opts, oper);
toc;

% Size of solution factor.
fprintf('\nSize outReg.Z: %d x %d\n', ...
    size(outReg.Z, 1), size(outReg.Z, 2));

% Check residuals.
res0 = norm(eqn.C * eqn.C');
res1 = mess_res2_norms(outReg.Z, 'riccati', ...
    eqn, opts, oper, opts.radi, []) / res0;
res2 = abs(eigs(@(x) eqn.A_' * (outReg.Z * ((outReg.Z' * x))) ...
    + (outReg.Z * ((outReg.Z' * (eqn.A_ * x)))) ...
    + eqn.C' * (eqn.C * x) - outReg.K' * (outReg.K * x), ...
    n, 1, 'LM')) / res0;

fprintf(['Residual computation -- RADI: %e | ' ...
    'mess_res2_norms: %e | eigs: %e \n'], ...
    outReg.res(end), res1, res2);

% Print convergence behavior.
if istest
    if min(outReg.res) >= opts.radi.res_tol
        error('MESS:TEST:accuracy', 'unexpectedly inaccurate result');
    end
else
    figure(1);
    semilogy(outReg.res,'LineWidth', 3);
    title('0 = C^T*C + A^T*X + X*A - X*B*B^T*X');
    xlabel('number of iteration steps');
    ylabel('normalized residual norm');
    pause(1);
end

fprintf('\n');


%% Solve filter Riccati equation.
fprintf('Solve filter Riccati equation\n');
fprintf('-----------------------------\n');

eqn.type     = 'N';
opts.LDL_T   = 1;
opts.radi.Z0 = ZC0;
opts.radi.Y0 = YC0;

% Some additional scaling for non-trivial LDL' term.
eqn.S        = diag(1:n_unstable);
eqn.B        = eqn.B * diag(1 ./ sqrt(1:n_unstable));

% Set corresponding residual.
opts.radi.W0 = WC0 * diag(1 ./ sqrt(1:n_unstable));
opts.radi.S0 = diag(1:n_unstable);

tic;
[outFil, eqn, opts, oper] = mess_lrradi(eqn, opts, oper);
toc;

% Size of solution factor.
fprintf('\nSize outFil.Z: %d x %d\n', ...
    size(outFil.Z, 1), size(outFil.Z, 2));

% Check residuals.
res0       = norm((eqn.B * eqn.S) * eqn.B');
eqn.S_diag = diag(eqn.S);
res3       = mess_res2_norms(outFil.Z, 'riccati', ...
    eqn, opts, oper, opts.radi, outFil.D) / res0;
ZD         = outFil.Z * outFil.D;
res4       = abs(eigs(@(x) eqn.A_ * (ZD * ((outFil.Z' * x))) ...
    + (ZD * ((outFil.Z' * (eqn.A_' * x)))) ...
    + eqn.B * (eqn.S * (eqn.B' * x)) - (outFil.K)' * ((outFil.K) * x), ...
    n, 1, 'LM')) / res0;

fprintf(['Residual computations -- RADI: %e | ' ...
    'mess_res2_norms: %e | eigs: %e \n'], ...
    outFil.res(end), res3, res4);

% Print convergence behavior.
if istest
    if min(outFil.res) >= opts.radi.res_tol
        error('MESS:TEST:accuracy', 'unexpectedly inaccurate result');
    end
else
    figure(2);
    semilogy(outFil.res,'LineWidth', 3);
    title('0 = B*S*B^T + A*Y + Y*A^T - Y*C^T*C*Y');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
