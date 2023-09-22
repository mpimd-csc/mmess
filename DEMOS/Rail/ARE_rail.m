function ARE_rail(k, shifts, inexact, Galerkin, type, istest)
% Computes the optimal feedback via the low-rank Newton-ADI [1] and RADI
% [2] methods for the selective cooling of Steel profiles application
% described in [3,4,5].
%
% Usage: ARE_Rail(k,shifts,inexact,Galerkin,istest)
%
% Inputs:
%
% k           refinement level of the model to use
%             (0 - 5, i.e. 109 - 79841 Dofs)
%             (optional, defaults to 2, i.e. 1357 Dofs)
%
% shifts      ADI shift selection strategy. Possible values:
%              'heur'        Penzl's heuristic shifts
%              'Wachspress'  Wachspress shifts, optimally solving the dense
%                          shift selection problem.
%             (optional, defaults to 'heur')
%
% inexact     use inexact Newton method
%             (optional, defaults to 0, i.e. false)
%
% Galerkin    activate Galerkin projection acceleration in Newton method.
%             This supersedes inexact Newton selection, i.e, disables it in
%             case both are on.
%             (optional, defaults to 0, i.e. no Galerkin acceleration)
%
% type        selector for the type of equation solved:
%               'LQR'  classic LQR ARE with Q=I, R=I and opts.LDL_T=false
%               'Hinf' robust control type ARE with Q=I,
%                      R=diag([1,1,1,1,-4,-4,-4]) and opts.LDL_T=true
%               'BR'   ARE with Q=I and R=-0.2401 I and opts.LDL_T=true as
%                      appearing in bounded real balanced truncation
%             (optional, defaults to 'LQR')
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% References:
% [1] P. Benner, J.-R. Li, T. Penzl, Numerical solution of large-scale
%     Lyapunov equations, Riccati equations, and linear-quadratic optimal
%     control problems, Numer. Lin. Alg. Appl. 15 (9) (2008) 755–777.
%     https://doi.org/10.1002/nla.622
%
% [2] P. Benner, Z. Bujanović, P. Kürschner, J. Saak, RADI: A low-rank
%     ADI-type algorithm for large scale algebraic Riccati equations,
%     Numer. Math. 138 (2) (2018) 301–330.
%     https://doi.org/10.1007/s00211-017-0907-5
%
% [3] J. Saak, Effiziente numerische Lösung eines
%     Optimalsteuerungsproblems für die Abkühlung von Stahlprofilen,
%     Diplomarbeit, Fachbereich 3/Mathematik und Informatik, Universität
%     Bremen, D-28334 Bremen (Sep. 2003).
%     https://doi.org/10.5281/zenodo.1187040
%
% [4] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lecture Notes in Computational Science and Engineering,
%     Springer-Verlag, Berlin/Heidelberg, Germany, 2005, pp. 353–356.
%     https://doi.org/10.1007/3-540-27909-1_19
%
% [5] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
narginchk(0, 6);
if nargin < 1
    k = 2;
end
if nargin < 2
    shifts = 'heur';
end
if nargin < 3
    inexact = false;
end
if nargin < 4
    Galerkin = false;
end
if nargin < 5
    istest = false;
end
if nargin < 6
    type = 'LQR';
end

%% set operation
[oper, opts] = operatormanager(struct(), 'default');
%% Problem data

eqn = mess_get_linear_rail(k);
switch type
    case 'LQR'
        opts.LDL_T = false;
        eqn.R = eye(size(eqn.B, 2)); % only used in here

    case 'Hinf'
        opts.LDL_T = true;
        eqn.Q = eye(size(eqn.C, 1));
        eqn.R = diag([1, 1, 1, 1, -4, -4, -4]);

    case 'BR'
        opts.LDL_T = true;
        eqn.Q = eye(size(eqn.C, 1));
        eqn.R = -0.2401 * eye(size(eqn.B, 2));

end

if opts.LDL_T
    mytitle = '0 = C^T Q C + A^T X E + E^T X A - E^T X B R^{-1} B^T X E';
else
    mytitle = '0 = C^T C + A^T X E + E^T X A - E^T X BB^T X E';
end

%%
% First we run the Newton-ADI Method

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.res_tol = 1e-14;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;

eqn.type = 'T';

%%
% Heuristic shift parameters via basic Arnoldi
opts.shifts.num_desired = 25;
opts.shifts.num_Ritz = 50;
opts.shifts.num_hRitz = 25;
n = oper.size(eqn, opts);
opts.shifts.b0 = ones(n, 1);
switch lower(shifts)

    case 'heur'
        opts.shifts.method = 'heur';

    case 'wachspress'
        opts.shifts.method = 'wachspress';
        opts.shifts.wachspress = 'T';

    case 'projection'
        opts.shifts.method = 'projection';
end
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 20;
opts.nm.res_tol = 1e-12;
opts.nm.rel_diff_tol = 1e-16;
if istest
    opts.nm.info = 0;
else
    opts.nm.info = 1;
end

opts.nm.accumulateRes = true;
opts.norm = 'fro';

if Galerkin
    opts.nm.linesearch = false;
    opts.nm.inexact = false;
    opts.nm.projection.freq = 2;
    opts.nm.projection.ortho = true;
elseif inexact
    opts.nm.linesearch = true;
    opts.nm.inexact = 'quadratic';
    opts.nm.projection.freq = 0;
    opts.nm.projection.ortho = false;
else
    opts.nm.linesearch = false;
    opts.nm.inexact = false;
    opts.nm.projection.freq = 0;
    opts.nm.projection.ortho = false;
end
opts.nm.res = struct('maxiter', 10, 'tol', 1e-6, 'orth', 0);
%%
t_mess_lrnm = tic;
outnm = mess_lrnm(eqn, opts, oper);
t_elapsed1 = toc(t_mess_lrnm);
if not(istest)
    mess_fprintf(opts, 'mess_lrnm took %6.2f seconds \n', t_elapsed1);
end

if istest
    if min(outnm.res) >= opts.nm.res_tol
        mess_err(opts, 'TEST:accuracy', 'unexpectedly inaccurate result');
    end
else
    figure(1);
    semilogy(outnm.res, 'LineWidth', 3);
    title(mytitle);
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
    mess_fprintf(opts, 'size outnm.Z: %d x %d\n\n', ...
                 size(outnm.Z, 1), size(outnm.Z, 2));
end

%%
% Lets try the RADI method and compare

% RADI-MESS settings
opts.shifts.history = opts.shifts.num_desired * size(eqn.C, 1);
opts.shifts.num_desired = opts.shifts.num_desired;

% choose either of the three shift methods, here
opts.shifts.method = 'gen-ham-opti';
% opts.shifts.method = 'heur';
% opts.shifts.method = 'projection';

opts.shifts.naive_update_mode = false; % .. Suggest false
% (smart update is faster;
%  convergence is the same).
opts.radi.compute_sol_fac = true;
opts.radi.get_ZZt = true;
opts.radi.maxiter = opts.adi.maxiter;
opts.norm = 2;
opts.radi.res_tol = opts.nm.res_tol;
opts.radi.rel_diff_tol = 0;
if istest
    opts.radi.info = 0;
else
    opts.radi.info = 1;
end

t_mess_lrradi = tic;
outradi = mess_lrradi(eqn, opts, oper);
t_elapsed2 = toc(t_mess_lrradi);

if not(istest)
    mess_fprintf(opts, 'mess_lrradi took %6.2f seconds \n', t_elapsed2);
end

if istest
    if min(outnm.res) >= opts.nm.res_tol
        mess_err(opts, 'TEST:accuracy', 'unexpectedly inaccurate result');
    end
else
    figure(2);
    semilogy(outradi.res, 'LineWidth', 3);
    title(mytitle);
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    mess_fprintf(opts, 'size outradi.Z: %d x %d\n\n', ...
                 size(outradi.Z, 1), size(outradi.Z, 2));
end

%% compare
if not(istest)
    figure(3);
    ls_nm = [outnm.adi.niter];
    ls_radi = 1:outradi.niter;

    semilogy(cumsum(ls_nm), outnm.res, 'k--', ...
             ls_radi, outradi.res, 'b-', ...
             'LineWidth', 3);

    title(mytitle);
    xlabel('number of solves with A+p*M');
    ylabel('normalized residual norm');
    legend('LR-NM', 'RADI');
end

if opts.LDL_T
    [ZN, DN] = mess_column_compression(outnm.Z, 'N', outnm.D, eps, not(istest));
    ZN = ZN * diag(sqrt(diag(DN)));

    [ZR, DR] = mess_column_compression(outradi.Z, 'N', outradi.D, eps, ...
                                       not(istest));
    ZR = ZR * diag(sqrt(diag(DR)));

else
    ZN = outnm.Z;
    ZR = outradi.Z;
end

D = blkdiag(eye(size(ZN, 2)), -eye(size(ZR, 2)));

Z = [ZN ZR];
f = @(x) (Z * (D * (Z' * x)));

if k < 4
    if exist('icare', 'file')
        X = icare(full(eqn.A_), eqn.B, eqn.C' * eqn.C, eqn.R, [], full(eqn.E_), []);
    elseif exist('care', 'file')
        X = care(full(eqn.A_), eqn.B, eqn.C' * eqn.C, eqn.R, [], full(eqn.E_));
    else
        X = ZN * ZN';
    end

    relerr = norm(ZN * ZN' - ZR * ZR') / norm(ZN * ZN');
    relerr(2) = max(abs(eigs(f, size(Z, 1)))) / norm(ZN' * ZN);
    relerr(3) = norm(X - ZN * ZN') / norm(X);
    relerr(4) = norm(X - ZR * ZR') / norm(X);
else

    relerr = max(abs(eigs(f, size(Z, 1)))) / norm(ZN' * ZN);

end

if istest
    if any(relerr(1:2) > 1e-6)
        % for large examples relerr computations appear to become unstable,
        % hence the comparably large tolerance
        mess_fprintf(opts, '%s %g %g %g %g\n', type, relerr);
        mess_err(opts, 'TEST:accuracy', ...
                 'Newton and RADI solution approximations deviate too much.');

    end
else
    if k < 4
        mess_fprintf(opts, ...
                     ['relative deviation of the two solution ', ...
                      'approximations:\n%g (dense 2 norm) ', ...
                      '%g (low-rank approx 2-norm)\n'], ...
                     relerr(1), relerr(2));
        if exist('icare', 'file') || exist('care', 'file')
            mess_fprintf(opts, ...
                         ['relative deviation of the two solution ', ...
                          'approximations from (i)care''s solution:\n', ...
                          '%g (NM) %g (RADI)\n'], ...
                         relerr(3), relerr(4));
        end
    else
        mess_fprintf(opts, ['relative deviation of the two solution ', ...
                            'approximations: % g\n'], ...
                     relerr);
    end
end
