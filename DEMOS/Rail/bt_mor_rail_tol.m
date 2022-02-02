function [Ar, Br, Cr] = bt_mor_rail_tol(k,tol,shifts,istest)
% bt_mor_FDM_tol computes a reduced order model via the standard Lyapunov
% balanced truncation (see e.g. [1]) for a finite difference discretized
% convection diffusion model on the unit square described in [2].
%
% Usage:
%    [Ar, Br, Cr] = bt_mor_rail_tol(k,tol,max_ord,n0,test)
%
% Inputs
%
% k           refinement level of the model to use
%             (0 - 5, i.e. 109 - 79841 Dofs)
%             (optional, defaults to 3, i.e. 5177 Dofs)
%
% tol         truncation tolerance for the Hankel singular values
%             (optional; defaults to 1e-6)
%
% shifts      shift selection used in ADI;  possible choices:
%               'heur'       :   Penzl heuristic shifts
%               'projection' :   projection shifts using the last columns
%                                of the solution factor
%             (optional, defaults to 'heur')
%
% istest      flag to determine whether this demo runs as a CI test or
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% Outputs
%
% Ar, Br, Cr  the reduced order system matrices.
%
% References
% [1] A. C. Antoulas, Approximation of Large-Scale Dynamical Systems, Vol.
%     6 of Adv. Des. Control, SIAM Publications, Philadelphia, PA, 2005.
%     https://doi.org/10.1137/1.9780898718713
%
% [2] J. Saak, Effiziente numerische Lösung eines
%     Optimalsteuerungsproblems für die Abkühlung von Stahlprofilen,
%     Diplomarbeit, Fachbereich 3/Mathematik und Informatik, Universität
%     Bremen, D-28334 Bremen (Sep. 2003).
%     https://doi.org/10.5281/zenodo.1187040
%
% [3] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353–356. https://doi.org/10.1007/3-540-27909-1_19
%
% [4] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%
%%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

narginchk(0,4);

% BT tolerance and maximum order for the ROM
if nargin<1, k=3; end
if nargin<2, tol=1e-5; end
if nargin<3, shifts='heur'; end
if nargin<4, istest=0; end


% ADI tolerance and maximum iteration number
opts.adi.maxiter = 100;              % maximum iteration number
opts.adi.res_tol = 1e-10;            % residual norm tolerance
opts.adi.rel_diff_tol = 1e-16;       % relative change norm tolerance
opts.adi.info = 1;                   % turn output on
opts.norm = 'fro';                   % Frobenius norm for stopping criteria

oper = operatormanager('default');
%%  Problem data
eqn = mess_get_linear_rail(k);              % load system matrices
n = oper.size(eqn, opts);                   % number of equations
%% Shift Parameters
opts.shifts.num_desired=25;                 % number of parameters for
                                            % 'heur' and 'wachspress'
switch lower(shifts)

    case 'heur'
        opts.shifts.method = 'heur';
        opts.shifts.num_Ritz=50;            % number Arnoldi steps with F
        opts.shifts.num_hRitz=25;           % Arnoldi steps with inv(F)
        opts.shifts.b0=ones(n,1);           % initial guess for Arnoldi

    case 'wachspress'
        opts.shifts.method = 'wachspress';
        opts.shifts.wachspress = 'T';

    case 'projection'
        opts.shifts.method = 'projection';
end
%% Compute low-rank factor of Controllability Gramian
eqn.type='N';                       % Lyapunov eq. for Controllability Gram.
t_mess_lradi = tic;
outB = mess_lradi(eqn, opts, oper); % run ADI iteration
t_elapsed1 = toc(t_mess_lradi);
fprintf(1,'mess_lradi took %6.2f seconds \n',t_elapsed1);

% residual norm plot
if istest
    if min(outB.res)>=opts.adi.res_tol
        error('MESS:TEST:accuracy','unexpectedly inaccurate result');
    end
else
    figure(1);
    semilogy(outB.res,'LineWidth',3);
    title('A X E^T + E X A^T = -BB^T');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outB.Z:');
disp(size(outB.Z));

%% Compute low-rank factor of Observability Gramian
eqn.type = 'T';                     % Lyapunov eq. for Observability Gram.
t_mess_lradi = tic;
outC = mess_lradi(eqn, opts, oper); % run ADI iteration
t_elapsed2 = toc(t_mess_lradi);
fprintf(1,'mess_lradi took %6.2f seconds \n',t_elapsed2);

% residual norm plot
if istest
    if min(outC.res)>=opts.adi.res_tol
        error('MESS:TEST:accuracy','unexpectedly inaccurate result');
    end
else
    figure(2);
    semilogy(outC.res,'LineWidth',3);
    title('A^T X E + E^T X A = -C^T C');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outC.Z:');
disp(size(outC.Z));

%% Compute reduced system matrices
% Perform Square Root Method
opts.srm.tol=tol;
opts.srm.max_ord=n;
opts.srm.info=2;
[TL,TR,HSV] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);
%% compute ROM matrices
Ar = TL'*oper.mul_A(eqn, opts, 'N', TR, 'N');
Br = TL'*eqn.B;
Cr = eqn.C*TR;
Er = eye(size(Ar,1));
%% Plots
ROM.A=Ar;
ROM.E=Er;
ROM.B=Br;
ROM.C=Cr;
if istest
    opts.sigma.info=0;
else
    opts.sigma.info=2;
end

opts.sigma.fmin=-6;
opts.sigma.fmax=4;

out = mess_sigma_plot(eqn, opts, oper, ROM); err = out.err;

if istest
    if max(err)>tol
        error('MESS:TEST:accuracy','unexpectedly inaccurate result');
    end
else
    figure;
    semilogy(HSV,'LineWidth',3);
    title('Computed Hankel singular values');
    xlabel('index');
    ylabel('magnitude');
end
