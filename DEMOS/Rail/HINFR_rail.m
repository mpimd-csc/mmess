function HINFR_rail(k, istest)
% Computes the a robust suboptimal Hinf feedback via the low-rank
% Riccati iteration [1] using the RADI method [2] for the selective cooling of
% Steel profiles application described in [3,4,5].
% 
% Usage: HINFR_rail(k, istest)
%
% Inputs: 
% 
% k           refinement level of the model to use 
%             (1-4, i.e. 1357-79841Dofs)
%             (optinal, defaults to 1)
%
% istest      flag to determine whether this demo runs as a CI test or 
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% References:
% [1] A. Lanzon, Y. Feng, B. D. O. Anderson, An iterative algorithm to 
%     solve algebraic Riccati equations with an indefinite quadratic term,
%     2007 European Control Conference (ECC), pp. 3033--3039, 2007.
%     https://doi.org/10.23919/ecc.2007.7068239
%
% [2] P. Benner, Z. Bujanović, P. Kürschner, J. Saak, RADI: A low-rank
%     ADI-type algorithm for large scale algebraic Riccati equations, 
%     Numer. Math. 138 (2) (2018) 301–330. 
%     https://doi.org/10.1007/s00211-017-0907-5.
%
% [3] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).  
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
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
%               2009-2020
%

%%
narginchk(0,2);
if nargin<1, k=1; end
if nargin<2, istest=0; end

%% set operation
oper = operatormanager('default');

%% Problem data
eqn = getrail(k);
% Reformulate as normalized Hinf control problem.
eqn.B1 = eqn.B;
eqn.B2 = eqn.B;
eqn.C1 = eqn.C;
eqn = rmfield(eqn,'B');
eqn = rmfield(eqn,'C');


%% Optional parameters.
opts.norm                = 2;
% Shift selection settings.
opts.shifts.history = 42;
opts.shifts.num_desired      = 5;
% choose either of the three shift methods, here
opts.shifts.method = 'gen-ham-opti';
%opts.shifts.method = 'heur';
%opts.shifts.method = 'projection';

% RADI settings
opts.shifts.naive_update_mode = false; % .. Suggest false (smart update is faster; convergence is the same).
opts.radi.compute_sol_fac     = 1;
opts.radi.get_ZZt             = 1;
opts.radi.maxiter             = 200;
opts.radi.res_tol             = 1.0e-10;
opts.radi.rel_diff_tol        = 1.0e-16;
opts.radi.info                = 1;

% Riccati iteration settings.
opts.ri.riccati_solver = 'radi';
opts.ri.maxiter        = 10;
opts.ri.res_tol        = 1.0e-10;
opts.ri.rel_diff_tol   = 1.0e-16;
opts.ri.compres_tol    = 1.0e-16;
opts.ri.info           = 1;

%% Solve the equation.
eqn.type = 'T';
gam      = 10;
eqn.B1   = 1/gam * eqn.B1; 
tic;
out = mess_lrri(eqn, opts, oper);
toc;

%% Residual behavior.
if istest
    if min(out.res) >= opts.ri.res_tol
       error('MESS:TEST:accuracy','unexpectedly inaccurate result'); 
   end
else
    figure(1);
    semilogy(out.res);
    hold on;
    for i = 1:length(out.radi), semilogy(out.radi(i).res); end
    hold off;
    title('0= C_1^TC_1 + A^TXM + M^TXA  + M^TX(\gamma^{-2}B_1B_1^T - B_2B_2^T)XM');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    legend('Riccati Iteration', 'RADI (step 1)', 'RADI (step 2)');
end
disp('size out.Z:');
disp(size(out.Z));
