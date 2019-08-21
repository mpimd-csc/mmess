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
%             (1-4, i.e. 1357-79841Dofs)
%             (optinal, defaults to 1)
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
%     https://doi.org/10.1002/nla.622.    
%
% [2] P. Benner, Z. Bujanović, P. Kürschner, J. Saak, RADI: A low-rank
%     ADI-type algorithm for large scale algebraic Riccati equations,
%     Numer. Math. 138 (2) (2018) 301–330.
%     https://doi.org/10.1007/s00211-017-0907-5.  
%
% [3] J. Saak, Effiziente numerische Lösung eines
%     Optimalsteuerungsproblems für die Abkühlung von Stahlprofilen,
%     Diplomarbeit, Fachbereich 3/Mathematik und Informatik, Universität
%     Bremen, D-28334 Bremen (Sep. 2003).   
%
% [4] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353–356. https://doi.org/10.1007/3-540-27909-1_19. 
%
% [5] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).  
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%


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

%%
narginchk(0,5);
if nargin<1, k=1; end
if nargin<2, shifts='heur'; end
if nargin<3, inexact=0; end
if nargin<4, Galerkin=0; end
if nargin<5, istest=0; end

%% set operation
oper = operatormanager('default');
%% Problem data

eqn=getrail(k);
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
tic;
outnm = mess_lrnm(eqn, opts, oper);
toc;

if istest
    if min(outnm.res)>=opts.nm.res_tol
       error('MESS:TEST:accuracy','unexpectedly innacurate result'); 
   end
else
    figure(1);
    disp(outnm.res);
    semilogy(outnm.res);
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

tic;
outradi = mess_lrradi(eqn, opts, oper);
toc;

if istest
    if min(outnm.res)>=opts.nm.res_tol
       error('MESS:TEST:accuracy','unexpectedly innacurate result'); 
   end
else
    figure(2);
    semilogy(outradi.res);
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

    semilogy(cumsum(ls_nm),outnm.res,'k--',ls_radi,outradi.res,'b-');
    title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
    xlabel('number of solves with A+p*M');
    ylabel('normalized residual norm');
    legend('LR-NM','RADI');
end