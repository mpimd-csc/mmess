function LQR_DAE3_SO(model,istest)
% Computes a Riccati feedback control for the constrained vibrational model
% from [1]
%
%
%   Usage:
%    LQR_DAE3_SO(model,test)
%
% Input:
%
%  model    choice of the example system. Possible values
%           'Stykel_small'   Stykel's mass spring damper system from [1]
%                            of dimension 600 in second order form (default)
%           'Stykel_large'   Stykel's mass spring damper system from [1]
%                            of dimension 6000 in second order form
%
%  istest    decides whether the function runs as an interactive demo or a
%            continuous integration test.
%            (optional; defaults to 0, i.e. interactive demo)
%
% References:
% [1] V. Mehrmann and T. Stykel. Balanced truncation model reduction for
%     large-scale systems in descriptor form.
%     In P. Benner, V. Mehrmann, and D. Sorensen, editors, Dimension
%     Reduction of Large-Scale Systems, volume 45 of Lecture Notes in
%     Computational Science and Engineering, pages 83–115. Springer-Verlag,
%     Berlin/Heidelberg, 2005. https://doi.org/10.1007/3-540-27909-1_3

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
narginchk(0,2);

if (nargin < 1)
    model = 'Stykel_small';
end

if ( nargin < 2 )
    istest = 0;
end
%% set operation
oper = operatormanager('dae_3_so');

%% load problem data
switch lower(model)
    case 'stykel_small'
             sys = load(sprintf('%s/../models/ms_ind3_by_t_stykel/g600.mat',...
                fileparts(mfilename('fullpath'))));
    case 'stykel_large'
            sys = load(sprintf('%s/../models/ms_ind3_by_t_stykel/g6000.mat',...
                fileparts(mfilename('fullpath'))));
end
eqn.M_=sys.M;
eqn.E_=sys.D;
eqn.K_=sys.K;
eqn.G_=sys.G;
eqn.haveE=1;
eqn.alpha = -0.02;
nv = size(eqn.M_,1);
np = size(eqn.G_,1);
eqn.B = sys.B(1:2*nv,:);
eqn.C = sys.C(:,1:2*nv);
eqn.type = 'T';

%%
% First we run the Newton-ADI Method
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

%% The actual Newton call
eqn.type = 'T';

t_mess_lrnm = tic;
outnm = mess_lrnm(eqn, opts, oper);
t_elapsed1 = toc(t_mess_lrnm);
fprintf(1,'mess_lrnm took %6.2f seconds \n',t_elapsed1);

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
fprintf('2-Norm of the resulting feedback matrix: %g\n', norm(outnm.Z,2));

disp('size outnm.Z:');
disp(size(outnm.Z));

%% Lets try the RADI metod and compare
% RADI-MESS settings
opts.shifts.history = opts.shifts.num_desired*size(eqn.C,1);
opts.shifts.num_desired=25;
opts.shifts.num_Ritz=50;
opts.shifts.num_hRitz=25;
opts.shifts.b0=ones(2*nv+np,1);

% choose either of the three shift methods, here
%opts.shifts.method = 'gen-ham-opti';
opts.shifts.method = 'heur';
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
fprintf(1,'mess_lrradi took %6.2f seconds \n' , t_elapsed2);

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

disp('size outradi.Z:');
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
