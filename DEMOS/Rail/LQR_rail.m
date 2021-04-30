function LQR_rail(k,shifts,inexact,Galerkin,istest)
% Computes the optimal feedback via the low-rank Newton-ADI [1] and RADI
% [2] methods for the selective cooling of Steel profiles application
% described in [3,4,5].
%
% Usage: LQR_Rail(k,shifts,inexact,Galerkin,istest)
%
% Inputs:
%
% k           refinement level of the model to use
%             (0 - 5, i.e. 109 - 79841 Dofs)
%             (optinal, defaults to 2, i.e. 1357 Dofs)
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
%             This superseeds inexact Newton selection, i.e, disables it in
%             case both are on.
%             (optional, defaults to 0, i.e. no Galerkin acceleration)
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
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353–356. https://doi.org/10.1007/3-540-27909-1_19
%
% [5] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%%
narginchk(0,5);
if nargin<1, k=2; end
if nargin<2, shifts='heur'; end
if nargin<3, inexact=0; end
if nargin<4, Galerkin=0; end
if nargin<5, istest=0; end

%% set operation
oper = operatormanager('default');
%% Problem data

eqn=mess_get_linear_rail(k);
%%
% First we run the Newton-ADI Method


% ADI tolerances and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.res_tol = 1e-14;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 1;

eqn.type = 'T';


%%
%Heuristic shift parameters via basic Arnoldi
opts.shifts.num_desired=25;
opts.shifts.num_Ritz=50;
opts.shifts.num_hRitz=25;
n=oper.size(eqn, opts);
opts.shifts.b0=ones(n,1);
switch lower(shifts)

    case 'heur'
        opts.shifts.method = 'heur';

    case 'wachspress'
        opts.shifts.method = 'wachspress';
        opts.shifts.wachspress = 'T';
end
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 8;
opts.nm.res_tol = 1e-10;
opts.nm.rel_diff_tol = 1e-16;
opts.nm.info = 1;
opts.nm.accumulateRes = 1;
opts.norm = 'fro';

if Galerkin
    opts.nm.linesearch = 0;
    opts.nm.inexact = 0;
    opts.nm.projection.freq=2;
    opts.nm.projection.ortho=1;
elseif inexact
    opts.nm.linesearch = 1;
    opts.nm.inexact = 'quadratic';
    opts.nm.projection.freq=0;
    opts.nm.projection.ortho=0;
else
    opts.nm.linesearch = 0;
    opts.nm.inexact = 0;
    opts.nm.projection.freq=0;
    opts.nm.projection.ortho=0;
end
opts.nm.res=struct('maxiter',10,'tol',1e-6,'orth',0);
%%
t_mess_lrnm = tic;
outnm = mess_lrnm(eqn, opts, oper);
t_elapsed1 = toc(t_mess_lrnm);
fprintf(1,'mess_lrnm took %6.2f seconds \n' , t_elapsed1);

if istest
    if min(outnm.res)>=opts.nm.res_tol
       error('MESS:TEST:accuracy','unexpectedly inaccurate result');
   end
else
    figure(1);
    disp(outnm.res);
    semilogy(outnm.res,'linewidth',3);
    title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end

disp('size outnm.Z:');
disp(size(outnm.Z));

%%
% Lets try the RADI method and compare

% RADI-MESS settings
opts.shifts.history = opts.shifts.num_desired*size(eqn.C,1);
opts.shifts.num_desired = opts.shifts.num_desired;

% choose either of the three shift methods, here
opts.shifts.method = 'gen-ham-opti';
%opts.shifts.method = 'heur';
%opts.shifts.method = 'projection';

opts.shifts.naive_update_mode = false; % .. Suggest false (smart update is faster; convergence is the same).
opts.radi.compute_sol_fac = 1;
opts.radi.get_ZZt = 1;
opts.radi.maxiter = opts.adi.maxiter;
opts.norm = 2;
opts.radi.res_tol = opts.nm.res_tol;
opts.radi.rel_diff_tol = 0;
opts.radi.info = 1;

t_mess_lrradi = tic;
outradi = mess_lrradi(eqn, opts, oper);
t_elapsed2 = toc(t_mess_lrradi);
fprintf(1,'mess_lrradi took %6.2f seconds \n', t_elapsed2);

if istest
    if min(outnm.res)>=opts.nm.res_tol
       error('MESS:TEST:accuracy','unexpectedly inaccurate result');
   end
else
    figure(2);
    semilogy(outradi.res,'linewidth',3);
    title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
end
disp('size outradi.Z:');
disp(size(outradi.Z));

%% compare
if not(istest)
    figure(3);
    ls_nm=[outnm.adi.niter];
    ls_radi=1:outradi.niter;

    semilogy(cumsum(ls_nm),outnm.res,'k--',ls_radi,outradi.res,'b-','linewidth',3);
    title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
    xlabel('number of solves with A+p*M');
    ylabel('normalized residual norm');
    legend('LR-NM','RADI');
end
