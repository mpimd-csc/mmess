function BT_sym_TripleChain(variant,istest)
%
% Computes a reduced order model (ROM) for the triple chain example of
% Truhar and Veselic [1] via Balanced truncation, e.g. [2], exploiting
% symmetry of the control system, i.e. equality of the system Gramians.
% In comparison to BT_TripleChain this demo exploits, that the first order
% form of the system has symmetric E and A and B=C', such that the two
% Lyapunov equations coincide.
%
% Usage:   BT_sym_TripleChain(version)
%
% Input:
%
% variant  Decides the Balanced Truncation version to use.
%          Possible values:
%          'FO' for reduction of the first order form to first order form
%          'VV' velocity-velocity balancing of the second order form to
%               second order form.
%          'PP' position-position balancing of the second order form to
%               second order form.
%          'PV' position-velocity balancing of the second order form to
%               second order form.
%          'VP' velocity-position balancing of the second order form to
%               second order form.%
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% References:
%
% [1] N. Truhar and K. Veselic, An efficient method for estimating the
%     optimal dampers’ viscosity for linear vibrating systems using
%     Lyapunov equation, SIAM J. Matrix Anal. Appl., 31 (2009), pp. 18–39.
%     https://doi.org/10.1137/070683052
%
% [2] A. C. Antoulas, Approximation of Large-Scale Dynamical Systems, Vol.
%     6 of Adv. Des. Control, SIAM Publications, Philadelphia, PA, 2005.
%     https://doi.org/10.1137/1.9780898718713

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%%
narginchk(0,2);

if nargin==0
    variant = 'FO';
end
if nargin<2
    istest=0;
end
format long e;
% set operation
oper = operatormanager('so_2');
% Problem data

n1=500;
alpha=.002;
Beta=alpha;
v=5;

[eqn.M_,eqn.E_,eqn.K_]=triplechain_MSD(n1,alpha,Beta,v);

s  = size(eqn.K_,1);

B = ones(s,1);
O = zeros(s,1);
Cp = ones(1,s);
Cv = O';
eqn.B = [B; O];
eqn.C = [Cp, Cv];

eqn.haveE=1;

%%
% ADI tolerances and maximum iteration number
opts.adi.maxiter = 300;
opts.adi.res_tol = 1e-10;
opts.norm = 'fro';

%opts.adi.rel_diff_tol = 1e-16;
opts.adi.rel_diff_tol = 0;
opts.adi.info = 1;
%opts.adi.accumulateK = 0;
%opts.adi.accumulateDeltaK = 0;
opts.adi.compute_sol_fac = 1;

eqn.type = 'N';
%%
%Heuristic shift parameters via projection
opts.shifts.num_desired=5;

opts.shifts.info=0;
opts.shifts.method = 'projection';

%%
% Compute Gramian Factor (one is enough, since the Gramians are equal)
t_mess_lradi = tic;
outB = mess_lradi(eqn, opts, oper);
t_elapsed = toc(t_mess_lradi);
fprintf(1,'mess_lradi took %6.2f seconds \n',t_elapsed);

if istest
    if min(outB.res)>=1e-1
       error('MESS:TEST:accuracy','unexpectedly inaccurate result');
   end
else
    figure(1);
    semilogy(outB.res,'LineWidth',3);
    title('0 = A X ^T + E X A^T - BB^T'); 
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end

disp('size outB.Z:');
disp(size(outB.Z));

%%
switch upper(variant)
    case 'FO'
        %%
        % Compute first order ROM
        opts.srm.max_ord = 150;
        opts.srm.tol = eps;
        opts.srm.info = 1;

        [TL,TR] = mess_square_root_method(eqn,opts,oper,outB.Z,outB.Z);

        ROM.E = eye(size(TL,2));
        ROM.A = TL'*oper.mul_A(eqn, opts, 'N', TR, 'N');
        ROM.B = TL'*eqn.B;
        ROM.C = eqn.C*TR;
        ROM.D = [];
    case 'VV'
        U = outB.Z(1:s,:);
        V = outB.Z(1:s,:);
    case 'PP'
        U = outB.Z(s+1:end,:);
        V = outB.Z(s+1:end,:);
    case 'PV'
        U = outB.Z(s+1:end,:);
        V = outB.Z(1:s,:);
    case 'VP'
        U = outB.Z(1:s,:);
        V = outB.Z(s+1:end,:);
end
if not(strcmp(variant,'FO'))
    max_ord = 75;
    tol = eps;
    inform = 1;

    [TL,TR] = square_root_method_SO(eqn.M_, max_ord, tol, inform, U, V);

    ROM.M = eye(size(TL,2));
    ROM.E = TL'*(eqn.E_*TR);
    ROM.K = TL'*(eqn.K_*TR);
    ROM.B = TL'*B;
    ROM.Cv = Cv*TR;
    ROM.Cp = Cp*TR;
end
%%
% plot results
opts.sigma.fmin = 1e-4;
opts.sigma.fmax = 1e0;
opts.sigma.nsample = 200;
if istest
    opts.sigma.info = 1;
else
    opts.sigma.info = 2;
end
out = mess_sigma_plot(eqn, opts, oper, ROM); err = out.err;
if istest
    if max(err) > 1000
        error('MESS:TEST:accuracy','unexpectedly inaccurate result %g', max(err));
    end
end

