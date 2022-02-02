function LQR_DAE2(problem, level, re, istest)
% Computes a stabilizing feedback by implicitly solving the generalized
% Riccati equation for the equivalent projected system on the hidden manifold.
%
% Inputs:
% problem       either 'Stokes' or 'NSE' to choose the Stokes demo or the
%               linearized Navier-Stokes-Equation.
%               (optional, defaults to 'Stokes')
%
% level         discretization level 1 through 5
%               (optional, only used in 'NSE' case, default: 1)
%
% re            Reynolds number 300, 400, or 500
%               (optional, only used in 'NSE' case, default: 500)
%
% istest        flag to determine whether this demo runs as a CI test or
%               interactive demo
%               (optional, defaults to 0, i.e. interactive demo)
%
% Note that the 'NSE' option requires additional data available in a
% separate 270MB archive and at least the 5th discretization level needs a
% considerable amount of main memory installed in your machine.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Set operations
oper = operatormanager('dae_2');

%% Problem data
if nargin<1, problem='stokes'; end
if nargin<2, level=1; end
if nargin<3, re=500; end
if nargin<4, istest=0; end

switch lower(problem)
    case 'stokes'
        nin = 5;
        nout = 5;
        nx = 10;
        ny = 10;
        [eqn.E_,eqn.A_,eqn.B,eqn.C]=stokes_ind2(nin,nout,nx,ny);
        eqn.haveE=1;
        st=trace(eqn.E_); % Stokes is FDM discretized, so so this is
        % the dimension of the velocity space
        eqn.st=st;
        eqn.B=eqn.B(1:st,:);
        eqn.C=eqn.C(:,1:st);
    case 'nse'
        [eqn, K0, ~] = mess_get_NSE(re, level);
        opts.nm.K0 = K0;
        opts.radi.K0 = K0;
    otherwise
        error('input ''problem'' must be either ''NSE'' or ''Stokes''');
end
%%
% First we run the Newton-ADI Method
opts.norm = 2;

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 300;
opts.adi.res_tol = 1e-12;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 1;
opts.adi.LDL_T=0;
eqn.type='T';
%%
n=size(eqn.A_, 1);
opts.shifts.num_desired=5;%*nout;
opts.shifts.num_Ritz=50;
opts.shifts.num_hRitz=25;
opts.shifts.method = 'projection';
opts.shifts.b0=ones(n,1);
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 20;
opts.nm.res_tol = 1e-10;
opts.nm.rel_diff_tol = 1e-16;
opts.nm.info = 1;
opts.nm.projection.freq=0;
opts.nm.projection.ortho=1;
%opts.nm.projection.meth='care_nwt_fac';
opts.nm.res=struct('maxiter',10,'tol',1e-6,'orth',0);
opts.nm.linesearch = 1;
opts.nm.inexact = 'superlinear';
opts.nm.tau = 0.1;
opts.nm.accumulateRes = 1;

%% use low-rank Newton-Kleinman-ADI
t_mess_lrnm = tic;
outnm = mess_lrnm(eqn, opts, oper);
t_elapsed1 =toc(t_mess_lrnm);
fprintf(1,'mess_lrnm took %6.2f seconds \n' , t_elapsed1);
if not(istest)
    figure(1);
    disp(outnm.res);
    semilogy(outnm.res,'LineWidth',3);
    title('0 = C^T C + A^T X E + E^T X A - E^T X BB^T X E');
    xlabel('number of newton iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outnm.Z:');
disp(size(outnm.Z));

%% Lets try the RADI method and compare
opts.norm = 2;
% RADI-MESS settings
opts.shifts.history = opts.shifts.num_desired*size(eqn.C,1);
opts.shifts.num_desired = opts.shifts.num_desired;

% choose either of the three shift methods, here
opts.shifts.method = 'gen-ham-opti';
%     opts.shifts.method = 'heur';
%     opts.shifts.method = 'projection';

opts.shifts.naive_update_mode = false; % .. Suggest false (smart update is faster; convergence is the same).
opts.shifts.info = 0;
opts.radi.compute_sol_fac = 1; % Turned on for numerical stability reasons.
opts.radi.get_ZZt = 0;
opts.radi.maxiter = opts.adi.maxiter;
opts.radi.res_tol = opts.nm.res_tol;
opts.radi.rel_diff_tol = 0;
opts.radi.info = 1;


t_mess_lrradi = tic;
outradi = mess_lrradi(eqn, opts, oper);
t_elapsed2 = toc(t_mess_lrradi);
fprintf(1,'mess_lrradi took %6.2f seconds \n' ,t_elapsed2);
if not(istest)
    figure();
    semilogy(outradi.res,'LineWidth',3);
    title('0 = C^TC + A^T X E + E^T X A - E^T X BB^T X E');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
end

%% compare
if istest
    if min(outnm.res)>=opts.nm.res_tol, error('MESS:TEST:accuracy','unexpectedly inaccurate result'); end
    if min(outradi.res)>=opts.radi.res_tol, error('MESS:TEST:accuracy','unexpectedly inaccurate result'); end
else
    figure();
    ls_nm=[outnm.adi.niter];
    ls_radi=1:outradi.niter;

    semilogy(cumsum(ls_nm),outnm.res,'k--',ls_radi,outradi.res,'b-','LineWidth',3);
    title('0 = C^T C + A^T X E + E^T X A - E^T X BB^T X E');
    xlabel('number of solves with A+p*M');
    ylabel('normalized residual norm');
    legend('LR-NM','RADI');
end
