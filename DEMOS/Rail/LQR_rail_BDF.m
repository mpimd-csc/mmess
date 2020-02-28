function LQR_rail_BDF(k)
% Computes the optimal feedback via the low-rank Rosenbrock[1,2] methods for
% the selective cooling of Steel profiles application described in [3,4,5].
% 
% Usage: LQR_Rail(k,shifts,inexact,Galerkin,istest)
%
% Inputs: 
% 
% k           k-step BDF method  
%             possible values: 1, ..., 6
%             (optinal, defaults to 2)
%
% References:
%
% [1] N. Lang, H. Mena, J. Saak, On the benefits of the LDLT factorization
%     for largescale differential matrix equation solvers, Linear Algebra
%     Appl. 480 (2015) 4471.  https://doi.org/10.1016/j.laa.2015.04.006.
%
% [2] N. Lang, Numerical methods for large-scale linear time-varying
%     control systems and related differential matrix equations,
%     Dissertation, Technische Universitt Chemnitz, Chemnitz, Germany,
%     logos-Verlag, Berlin, ISBN 978-3-8325-4700-4 (Jun. 2017).    
%     URL https://www.logos-verlag.de/cgi-bin/buch/isbn/4700
%
% [3] J. Saak, Effiziente numerische Lsung eines
%     Optimalsteuerungsproblems fr die Abkhlung von Stahlprofilen,
%     Diplomarbeit, Fachbereich 3/Mathematik und Informatik, Universitt
%     Bremen, D-28334 Bremen (Sep. 2003).   
%
% [4] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353356. https://doi.org/10.1007/3-540-27909-1_19. 
%
% [5] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universitt Chemnitz, Chemnitz, Germany (Jul. 2009).  
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
%               2009-2020
%
%%
narginchk(0,1);
if nargin<1, k =2; end
%%
% set operation
oper = operatormanager('default');
% Problem data

eqn=getrail(0);
%%
opts.norm = 'fro';

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.res_tol = 1e-14;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;
opts.adi.compute_sol_fac = 1;
opts.cc_info = 0;

eqn.type = 'T';


%%
%Heuristic shift parameters via basic Arnoldi 
n=oper.size(eqn, opts);
opts.shifts.num_desired=7;
opts.shifts.num_Ritz=50;
opts.shifts.num_hRitz=25;
opts.shifts.method = 'heur';

opts.shifts.b0=ones(n,1);
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 8;
opts.nm.res_tol = 1e-10;
opts.nm.rel_diff_tol = 1e-16;
opts.nm.info = 0;
opts.norm = 'fro';
opts.nm.accumulateRes = 1;
opts.nm.linesearch = 1;

%%
% BDF parameters
opts.bdf.time_steps = 0 : 50: 4500;
opts.bdf.step = k;
opts.bdf.info = 1;
opts.bdf.save_solution = 0;
opts.bdf.startup_iter = 7;
%%
tic;
[out_bdf]=mess_bdf_dre(eqn,opts,oper);
toc;
%%
y = zeros(1,length(out_bdf.Ks));
for i=1:length(out_bdf.Ks)
    y(i) = out_bdf.Ks{i}(1,77);
end
x = opts.bdf.time_steps;
figure(1);
plot(x,y);
title('evolution of component (1,77) of the optimal feedback');
