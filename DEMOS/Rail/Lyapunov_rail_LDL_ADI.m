function Lyapunov_rail_LDL_ADI(k,shifts,implicit,istest)
% Computes the solution of the generalized Lyapunov equation via the
% low-rank ADI iteration in both ZZ^T [1,2] and LDL^T[3] formulation for the
% selective cooling of Steel profiles application described in [4,5,6].
% 
% Usage:      Lyapunov_rail_LDL_ADI(k,shifts,istest)
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
%                            shift selection problem.
%              'projection'  projection shifts [7]
%             (optional, defaults to 'heur')
% 
% implicit    projection shifts with implicit computation of projected A
%             matrix
%
% istest      flag to determine whether this demo runs as a CI test or 
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% References:
%
% [1] P. Benner, J.-R. Li, T. Penzl, Numerical solution of large-scale
%     Lyapunov equations, Riccati equations, and linear-quadratic optimal
%     control problems, Numer. Lin. Alg. Appl. 15 (9) (2008) 755–777. 
%     https://doi.org/10.1002/nla.622.    
%
% [2] P. Kürschner, Efficient low-rank solution of large-scale matrix
%     equations, Dissertation, Otto-von-Guericke-Universität, Magdeburg,
%     Germany, shaker Verlag, ISBN 978-3-8440-4385-3 (Apr. 2016). 
%     URL http://hdl.handle.net/11858/00-001M-0000-0029-CE18-2
%
% [3] N. Lang, H. Mena, J. Saak, On the benefits of the LDLT factorization
%     for largescale differential matrix equation solvers, Linear Algebra
%     Appl. 480 (2015) 44–71.  
%     https://doi.org/10.1016/j.laa.2015.04.006.
%
% [4] J. Saak, Effiziente numerische Lösung eines
%     Optimalsteuerungsproblems für die Abkühlung von Stahlprofilen,
%     Diplomarbeit, Fachbereich 3/Mathematik und Informatik, Universität
%     Bremen, D-28334 Bremen (Sep. 2003).   
%
% [5] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353–356. https://doi.org/10.1007/3-540-27909-1_19. 
%
% [6] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).  
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%
% [7] P. Benner, P. Kürschner, J. Saak, Self-generating and efficient shift
%     parameters in ADI methods for large Lyapunov and Sylvester equations,
%     Electron. Trans. Numer. Anal. 43 (2014) 142–162.
%     URL http://etna.mcs.kent.edu/volumes/2011-2020/vol43/abstract.php?vol=43&pages=142-162

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
narginchk(0,4);
if nargin<1, k=1; end
if nargin<2, shifts='wachspress'; end
if nargin<3, implicit=0; end
if nargin<4, istest=0; end
%%
% set operation
oper = operatormanager('default');
% Problem data
eqn = getrail(k);
%%

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.res_tol = 1e-12;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 1;

eqn.type = 'T';


%%
%Heuristic shift parameters via basic Arnoldi 
n=oper.size(eqn, opts);
opts.shifts.num_Ritz=50;
opts.shifts.num_hRitz=25;
opts.shifts.num_desired = 6;

opts.shifts.b0=ones(n,1);
switch lower(shifts)
    case 'heur'
        opts.shifts.method = 'heur';
    case 'wachspress'
        opts.shifts.method = 'wachspress';
        opts.shifts.wachspress = 'T';
    case 'projection'
        opts.shifts.method = 'projection';
        opts.shifts.implicitVtAV = implicit;
end

opts.shifts.p=mess_para(eqn,opts,oper);
opts.norm = 'fro';

%%
fprintf('########################\n');
fprintf('# ADI with ZZ^T:        \n');
fprintf('########################\n');
tic;
out = mess_lradi(eqn, opts, oper);
toc;

if istest
    if min(out.res)>=opts.adi.res_tol
       error('MESS:TEST:accuracy','unexpectedly innacurate result'); 
   end
else
    figure(1);
    semilogy(out.res);
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size out.Z:');
disp(size(out.Z));

%% Set LDL fields
opts.LDL_T = 1;
eqn.S = diag([4,4,9,9,16,16]);
eqn.G = eqn.C' * diag([4,4,9,9,16,16].^(-0.5));

%%
fprintf('########################\n');
fprintf('# ADI with LDL^T:       \n');
fprintf('########################\n');
opts.shifts.p=mess_para(eqn,opts,oper);
tic;
out1 = mess_lradi(eqn, opts, oper);
toc;

if istest
    if min(out.res)>=opts.adi.res_tol
       error('MESS:TEST:accuracy','unexpectedly innacurate result'); 
   end
else
    figure(2);
    semilogy(out1.res);
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size out1.Z:');
disp(size(out1.Z));

%% Difference of Lyapunov solutions
if k<3
    err = norm(out.Z * out.Z' - out1.Z * out1.D * out1.Z') / norm(out.Z * out.Z');
    fprintf('Relative difference between solution with and without LDL^T: \t %g\n', err);
     if err>1e-14
         error('MESS:TEST:accuracy','unexpectedly innacurate result');
     end
end