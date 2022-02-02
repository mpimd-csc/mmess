function LQR_rail_Rosenbrock(k)
% Computes the optimal feedback via the low-rank Rosenbrock[1,2] methods for
% the selective cooling of Steel profiles application described in [3,4,5].
%
% Usage: LQR_Rail(k,shifts,inexact,Galerkin,istest)
%
% Inputs:
%
% k           k-stage Rosenbrock method
%             possible values: 1, 2
%             (optional, defaults to 2)
%
% References:
%
% [1] N. Lang, H. Mena, J. Saak, On the benefits of the LDLT factorization
%     for large-scale differential matrix equation solvers, Linear Algebra
%     Appl. 480 (2015) 44–71.  https://doi.org/10.1016/j.laa.2015.04.006
%
% [2] N. Lang, Numerical methods for large-scale linear time-varying
%     control systems and related differential matrix equations,
%     Dissertation, Technische Universität Chemnitz, Chemnitz, Germany,
%     logos-Verlag, Berlin, ISBN 978-3-8325-4700-4 (Jun. 2017).
%     URL https://www.logos-verlag.de/cgi-bin/buch/isbn/4700
%
% [3] J. Saak, Effiziente numerische Lösung eines
%     Optimalsteuerungsproblems für die Abkühlung von Stahlprofilen,
%     Diplomarbeit, Fachbereich 3/Mathematik und Informatik, Universität
%     Bremen, D-28334 Bremen (Sep. 2003).
%     https://doi.org/10.5281/zenodo.1187040
%
% [4] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353–356. https://doi.org/10.1007/3-540-27909-1_19
%
% [5] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

narginchk(0,1);
if nargin<1, k=2; end
%%
% set operation
oper = operatormanager('default');

% Problem data
eqn = mess_get_linear_rail(0);
%%
opts.norm = 'fro';
% ADI tolerances and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.res_tol = 1e-14;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;
opts.adi.compute_sol_fac = 1;
opts.adi.accumulateK = 1;

eqn.type = 'T';

%%
%Heuristic shift parameters via projection
opts.shifts.num_desired=7;
opts.shifts.method = 'projection';
%%
% Rosenbrock parameters
opts.rosenbrock.time_steps = 0 : 50 : 4500;
opts.rosenbrock.stage = k;
opts.rosenbrock.info = 1;
opts.rosenbrock.gamma = 1 + 1 / sqrt(2);
opts.rosenbrock.save_solution = 0;
%%
t_mess_rosenbrock_dre = tic;
[out_ros]=mess_rosenbrock_dre(eqn,opts,oper);
t_elapsed = toc(t_mess_rosenbrock_dre);
fprintf(1,'mess_rosenbrock_dre took %6.2f seconds \n',t_elapsed);

y = zeros(1,length(out_ros.Ks));
for i=1:length(out_ros.Ks)
    y(i) = out_ros.Ks{i}(1,77);
end
x = opts.rosenbrock.time_steps;
figure(1);
plot(x,y,'LineWidth',3);
title('evolution of component (1,77) of the optimal feedback');
xlabel('time');
ylabel('magnitude');
