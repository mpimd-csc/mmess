function bt_mor_DAE2(problem,lvl,re,istest)
% Computes a standard ROM by implicitly solving the generalized Lyapunov
% equations for the equivalent projected system on the hidden manifold.
%
% Inputs:
% problem       either 'Stokes' or 'NSE' to choose the Stokes demo or the
%               linearized Navier-Stokes-Equation.
%               (ooptional, defaults to 'Stokes')
%
% lvl           discretization level 1 through 5
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
% See:
% P. Benner, J. Saak, M. M. Uddin, Balancing based model reduction for
% structured index-2 unstable descriptor systems with application to flow
% control, Numerical Algebra, Control and Optimization 6 (1) (2016) 1–20.
% https://doi.org/10.3934/naco.2016.6.1

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% ADI tolerance and maximum iteration number
opts.adi.maxiter = 350;
opts.adi.res_tol = sqrt(eps);
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 1;
opts.shifts.info = 1;
opts.norm = 'fro';

oper = operatormanager('dae_2');
%% Problem data
if nargin<1, problem='stokes'; end
if nargin<2, lvl=1; end
if nargin<3, re=500; end
if nargin<4, istest=0; end

problem = lower(problem);

switch problem
    case 'stokes'
        nin = 5;
        nout = 5;
        nx = 10;
        ny = 10;
        [eqn.E_,eqn.A_,eqn.Borig,eqn.Corig]=stokes_ind2(nin,nout,nx,ny);
        n=size(eqn.E_,1);
        eqn.haveE=1;
        st=trace(eqn.E_); % Stokes is FDM discretized, so so this is
                          % the dimension of the velocity space
        eqn.st=st;
        eqn.B=eqn.Borig(1:st,:);
        eqn.C=eqn.Corig(:,1:st);
    case 'nse'
        [eqn, K0primal, K0dual] = mess_get_NSE( re, lvl);
        st=eqn.st;
        n=size(eqn.E_,1);
    otherwise
        error('input ''problem'' must be either ''NSE'' or ''Stokes''');
end
%%
eqn.type='N';
% Activate stabilizing (Bernoulli) feedback
if strcmp(problem,'nse')
    eqn.V=-K0primal';
    eqn.U = eqn.B;
    eqn.haveUV=1;
end


opts.shifts.num_desired=6;
opts.shifts.num_Ritz=40;
opts.shifts.num_hRitz=40;
opts.shifts.method='projection';

opts.shifts.b0=ones(size(eqn.A_,1),1);

t_mess_lradi = tic;
outB = mess_lradi(eqn,opts,oper);
t_elapsed1 = toc(t_mess_lradi);
fprintf(1,'mess_lradi took %6.2f seconds \n',t_elapsed1);

if not(istest)
    figure();
    semilogy(outB.res,'linewidth',3);
    title('AXM^T + MXA^T = -BB^T');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outB.Z:');
disp(size(outB.Z));

%%
eqn.type = 'T';
% Activate stabilizing (Bernoulli) feedback (for the dual system)
if strcmp(problem,'nse')
    eqn.U=-K0dual';
    eqn.V = eqn.C';
    eqn.haveUV=1;
end


opts.shifts.num_desired=6;
opts.shifts.num_Ritz=40;
opts.shifts.num_hRitz=40;
opts.shifts.method='projection';

opts.shifts.b0=ones(size(eqn.A_,1),1);

t_mess_lradi = tic;
outC = mess_lradi(eqn, opts, oper);
t_elapsed2 = toc(t_mess_lradi);
fprintf(1,'mess_lradi took %6.2f seconds \n',t_elapsed2);

if not(istest)
    figure();
    semilogy(outC.res,'linewidth',3);
    title('A^TXM + M^TXA = -C^TC');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outC.Z:');
disp(size(outC.Z));

%% Compute reduced system matrices
% Perform Square Root Method  (SRM)

% BT tolerance and maximum order for the ROM
t_SRM = tic;
opts.srm.tol=1e-5;
opts.srm.max_ord=250;

% SRM verbosity
if istest
    opts.srm.info=1;
else
    opts.srm.info=2;
end

%The actual SRM
[TL,TR,hsv] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);

%%
ROM.A = TL'*(eqn.A_(1:st,1:st)*TR);
ROM.B = TL'*eqn.B(1:st,:);
ROM.C = eqn.C(:,1:st)*TR;

t_elapsed3 = toc(t_SRM);
fprintf(1,'computation of reduced system matrices took %6.2f seconds \n',t_elapsed3);

%%
t_eval_ROM = tic;
%% Evaluate the ROM quality
% while the Gramians are computed exploiting the DAE structure, due to the
% construction of the function handles we can not do so for the transfer
% function. Therfore we need to extend the matrices B and C and call the
% 'default' usfs for unstructured computation:
switch lower(problem)
    case 'stokes'
        eqn.B=eqn.Borig;
        eqn.C=eqn.Corig;
    case 'nse'
        n = size(eqn.A_,1);
        eqn.B(st+1:n,:) = zeros(n-st,size(eqn.B,2));
        eqn.C(:,st+1:n) = zeros(size(eqn.C,1),n-st);
end
oper = operatormanager('default');

if istest
    opts.sigma.info=0;
else
    opts.sigma.info=2;
end

opts.sigma.fmin=-3;
opts.sigma.fmax=4;

out = mess_sigma_plot(eqn, opts, oper, ROM); err = out.err;

t_elapsed4 = toc(t_eval_ROM);
fprintf(1,'evaluation of rom quality took %6.2f seconds \n' ,t_elapsed4);

%%
if istest
    if max(err)>=opts.srm.tol, error('MESS:TEST:accuracy','unexpectedly inaccurate result'); end
else
    figure;
    semilogy(hsv,'linewidth',3);
    title('Computed Hankel singular values');
    xlabel('index');
    ylabel('magnitude');
end
%%
fprintf(['\nComputing open loop step response of original and reduced order ' ...
    'systems and time domain MOR errors\n']);
open_step(eqn,ROM.A,ROM.B,ROM.C,problem,istest);
%%
fprintf('\nComputing ROM based feedback\n');
if exist('care', 'file')
    [~,~,Kr]=care(ROM.A,ROM.B,ROM.C'*ROM.C,eye(size(ROM.B,2)));
else
    Y = care_nwt_fac([],ROM.A,ROM.B,ROM.C,1e-12,50);
    Kr = (Y*ROM.B)'*Y;
end
K=[Kr*TL'*eqn.E_(1:st,1:st),zeros(size(Kr,1),n-st)];
%%
fprintf(['\nComputing closed loop step response of original and reduced order ' ...
    'systems and time domain MOR errors\n']);
closed_step(eqn,ROM.A,ROM.B,ROM.C,problem,K,Kr,istest);

