function LQR_DAE2_SO(istest)
% Computes a Riccati feedback control for a constrained variant of the
% triple chain oscillator model from [1] via the Newton-ADI and RADI
% methods. The code is using the ideas for second order matrix exploitation
% described in [2,3] and second order DAEs from [4].
%
% Input:
% istest  decides whether the function runs as an interactive demo or a
%         continuous integration test. (optional; defaults to 0, i.e.
%         interactive demo)
%
% References
% [1] N. Truhar, K. Veselić, An efficient method for estimating the
%     optimal dampers’ viscosity for linear vibrating systems using
%     Lyapunov equation, SIAM J. Matrix Anal. Appl. 31 (1) (2009) 18–39.
%     https://doi.org/10.1137/070683052
%
% [2] P. Benner, J. Saak, Efficient Balancing based MOR for Second Order
%     Systems Arising in Control of Machine Tools, in: I. Troch,
%     F. Breitenecker (Eds.), Proceedings of the MathMod 2009, no. 35 in
%     ARGESIM-Reports, Vienna Univ. of Technology, ARGE Simulation News,
%     Vienna, Austria, 2009, pp. 1232–1243, https://doi.org/10.11128/arep.35
%
% [3] P. Benner, P. Kürschner, J. Saak, Improved second-order balanced
%     truncation for symmetric systems, IFAC Proceedings Volumes (7th
%     Vienna International Conference on Mathematical Modelling) 45 (2)
%     (2012) 758–762. https://doi.org/10.3182/20120215-3-AT-3016.00134
%
% [4] M. M. Uddin, Computational methods for model reduction of large-scale
%     sparse structured descriptor systems, Dissertation,
%     Otto-von-Guericke-Universität, Magdeburg, Germany (2015).
%     URL http://nbn-resolving.de/urn:nbn:de:gbv:ma9:1-6535

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%%
if nargin<1, istest=0; end

%% set operation
oper = operatormanager('dae_2_so');

%% load problem data
% generate problem
nv = 500; np=10; nin=2; nout=3; p=0.2;
alpha = 0.1; beta = 0.1; v = 5;

[M, D, K]=triplechain_MSD(nv,alpha,beta,v);
nv=size(M,1);
G = zeros(nv,np);
while rank(full(G))~=np, G = sprand(np,nv,p);end
eqn.M_=M;
eqn.E_=-D;
eqn.K_=-K;
eqn.G_=G;
eqn.haveE=1;
eqn.alpha = -0.02;
eqn.B = [zeros(nv, nin);rand(nv,nin)];
eqn.C = [zeros(nout,nv),rand(nout,nv)];
eqn.type = 'N';

%% condition numbers of generated input data
fprintf('Condition numbers of the generated input data\n');
As = full([zeros(nv,nv),eye(nv,nv),zeros(nv,np); ...
     K, D, G';...
     G,zeros(np,nv), zeros(np,np)]);

fprintf('cond(M)=%e\n',condest(eqn.M_));
fprintf('cond(D)=%e\n',condest(eqn.E_));
fprintf('cond(K)=%e\n',condest(eqn.K_));
fprintf('cond(A)=%e\n\n',condest(As));

%% options
opts.norm = 'fro';

% ADI options
opts.adi.maxiter = 1000;
opts.adi.res_tol = 1e-15;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;
opts.shifts.num_desired=25;
opts.shifts.num_Ritz=50;
opts.shifts.num_hRitz=25;
opts.shifts.b0=ones(2*nv+np,1);
opts.shifts.method = 'projection';
opts.shifts.num_desired=6;

% Newton options and maximum iteration number
opts.nm.maxiter = 20;
opts.nm.res_tol = 1e-10;
opts.nm.rel_diff_tol = 1e-16;
opts.nm.info = 1;
opts.nm.accumulateRes = 0;
opts.nm.linesearch = 1;
opts.nm.projection.freq=0;
opts.nm.projection.ortho=1;
opts.nm.res=struct('maxiter',10,'tol',1e-6,'orth',0);

%%
% the actual Newton call
eqn.type = 'T';
t_mess_lrnm = tic;
outnm = mess_lrnm(eqn, opts, oper);
t_elapsed1 = toc(t_mess_lrnm);
fprintf(1,'mess_lrnm took %6.2f seconds ',t_elapsed1);
if istest
    if min(outnm.res)>=opts.nm.res_tol
        error('MESS:TEST:accuracy','unexpectedly inaccurate result');
    end
else
    figure(1);
    semilogy(outnm.res,'linewidth',3);
    title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end

disp('outnm.Z:');
disp(size(outnm.Z));

%% Lets try the RADI method and compare
opts.norm = 2;
% RADI-MESS settings
opts.shifts.history = opts.shifts.num_desired*size(eqn.C,1);
opts.shifts.num_desired = 5;

% choose either of the three shift methods, here
%opts.shifts.method = 'gen-ham-opti';
%opts.shifts.method = 'heur';
opts.shifts.method = 'projection';

opts.shifts.naive_update_mode = false;
% .. Suggest false (smart update is faster; convergence is the same).
opts.radi.compute_sol_fac = 1;
opts.radi.maxiter = opts.nm.maxiter * opts.adi.maxiter;
opts.radi.get_ZZt = 1;
opts.radi.res_tol = opts.nm.res_tol;
opts.radi.rel_diff_tol = 0;
opts.radi.info = 1;

t_mess_lrradi = tic;
outradi = mess_lrradi(eqn, opts, oper);
t_elapsed2 = toc(t_mess_lrradi);
fprintf(1,'mess_lrradi took %6.2f seconds ' ,t_elapsed2);
if istest
    if min(outradi.res)>=opts.radi.res_tol
        error('MESS:TEST:accuracy','unexpectedly inaccurate result');
    end
else
    figure(2);
    semilogy(outradi.res,'linewidth',3);
    title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
end

disp('outradi.Z:');
disp(size(outradi.Z));

%% compare
if not(istest)
    figure(3);
    ls_nm=cumsum([outnm.adi.niter]);
    ls_radi=1:outradi.niter;

    semilogy(ls_nm,outnm.res,'k--',ls_radi,outradi.res,'b-','linewidth',3);
    title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
    xlabel('number of solves with A+p*M');
    ylabel('normalized residual norm');
    legend('LR-NM','RADI');
end
