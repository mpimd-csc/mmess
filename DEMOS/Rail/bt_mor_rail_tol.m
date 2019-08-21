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
%             (1-4, i.e. 1357-79841Dofs)
%             (optinal, defaults to 1)
%
% tol         truncation tolerance for the Hankel singular values
%             (optional; defalts to 1e-6)
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
% Ar, Br, Cr  the reduced orde system matrices.
%
% References
% [1] A. C. Antoulas, Approximation of Large-Scale Dynamical Systems, Vol.
%     6 of Adv. Des. Control, SIAM Publications, Philadelphia, PA, 2005.
%     https://doi.org/10.1137/1.9780898718713.  
%
% [2] J. Saak, Effiziente numerische Lösung eines
%     Optimalsteuerungsproblems für die Abkühlung von Stahlprofilen,
%     Diplomarbeit, Fachbereich 3/Mathematik und Informatik, Universität
%     Bremen, D-28334 Bremen (Sep. 2003).   
%
% [3] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353–356. https://doi.org/10.1007/3-540-27909-1_19. 
%
% [4] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).  
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%
%%

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
%               2009-2019
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
eqn = getrail(k);                           % load system matrices
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

opts.shifts.p=mess_para(eqn,opts,oper); % compute shift parameters

disp('(Initial) ADI shifts');
disp(opts.shifts.p);
%% Compute low-rank factor of Controllability Gramian
eqn.type='N';                       % Lyapunov eq. for Controllability Gram.
tic;
outB = mess_lradi(eqn, opts, oper); % run ADI iteration
toc;
% residual norm plot
if istest
    if min(outB.res)>=opts.adi.res_tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result');
    end
else
    figure(1);
    semilogy(outB.res);
    title('AXM^T + MXA^T = -BB^T');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outB.Z:');
disp(size(outB.Z));

%% Compute low-rank factor of Observatibility Gramian
eqn.type = 'T';                     % Lyapunov eq. for Observability Gram.
tic;
outC = mess_lradi(eqn, opts, oper); % run ADI iteration
toc;
% residual norm plot
if istest
    if min(outC.res)>=opts.adi.res_tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result');
    end
else
    figure(2);
    semilogy(outC.res);
    title('A^TXM + M^TXA = -C^TC');
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
[TL,TR,hsv] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);
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

opts.sigma.fmin=-3;
opts.sigma.fmax=4;

err = mess_sigma_plot(eqn, opts, oper, ROM);

if istest
    if max(err)>tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result');
    end
else
    figure;
    semilogy(hsv);
    title('Computed Hankel singular values');
    xlabel('index');
    ylabel('magnitude');
end