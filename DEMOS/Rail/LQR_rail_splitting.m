function out =  LQR_rail_splitting(k, exp_action, method,istest)

% Computes the optimal feedback via low-rank splitting schemes [1, 2] for
% the selective cooling of Steel profiles application described in [3,4,5].
%
% Inputs:
%
% k           refinement level of the model to use
%             (0 - 5, i.e. 109 - 79841 Dofs)
%             (optional, defaults to 2, i.e. 1357 Dofs)
%
% exp_action  parameters for computing matrix-exponential actions, see
%             mat-eqn-solvers/private/mess_exp_action.m for details
%             (optional, defaults to exp_action.method = 'Krylov',
%             exp_action.tol = 1e-8)
%
% method      choice of splitting scheme; structure with fields 'order',
%             'additive' and 'symmetric'
%             (optional, defaults to order = 2, additive = false, symmetric
%             = false)
%             NOTE: the additive schemes typically need a running parallel
%             pool in order to be competitive
%
%
% References:
%
% [1] T. Stillfjord, Low-rank second-order splitting of large-scale
%     differential Riccati equations, IEEE Trans. Autom. Control, 60
%     (2015), pp. 2791-2796. https://doi.org/10.1109/TAC.2015.2398889
%
% [2] T. Stillfjord, Adaptive high-order splitting schemes for large-scale
%     differential Riccati equations, Numer. Algorithms,  (2017).
%     https://doi.org/10.1007/s11075-017-0416-8
%
% [3] J. Saak, Effiziente numerische Lösung eines
%     Optimalsteuerungsproblems für die Abkühlung von
%     Stahlprofilen, Diplomarbeit, Fachbereich 3/Mathematik und Informatik,
%     Universität Bremen, D-28334 Bremen (Sep. 2003).
%     https://doi.org/10.5281/zenodo.1187040
%
% [4] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353-356. https://doi.org/10.1007/3-540-27909-1_19
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


if nargin < 1
    k = 2;
end
if nargin < 2
    exp_action.method = 'Krylov';
    exp_action.tol = 1e-8;
end
if nargin < 3
    method.order = 2;
    method.additive = false;
    method.symmetric = false;
end
if nargin < 4
    istest = 0;
end
%% Equation parameters
% Default (E, A, B, C) system
oper = operatormanager('default');

eqn = mess_get_linear_rail(k);
eqn.Rinv = 1;

eqn.type = 'T';

eqn.L0 = rand(size(eqn.A_, 1), 1);
eqn.D0 = 1;


%% General splitting parameters
opts.splitting.time_steps = 0 : 50 : 4500;
opts.splitting.order = method.order;
opts.splitting.additive = method.additive;
opts.splitting.symmetric = method.symmetric;
opts.splitting.info = 2;
opts.splitting.intermediates = 1;
opts.splitting.trunc_tol = 1e-10;

% Quadrature (for integral term) parameters
opts.splitting.quadrature.type = 'adaptive';
opts.splitting.quadrature.tol = 1e-4 ;

%% Matrix exponential action parameters
opts.exp_action = exp_action;

%% Compute the approximation
t_mess_splitting_dre = tic;
[out, ~, opts, ~] = mess_splitting_dre(eqn, opts, oper);
t_elapsed = toc(t_mess_splitting_dre);
fprintf(1,'mess_splitting_dre took %6.2f seconds \n',t_elapsed);

%%
if not(istest)
    t = opts.splitting.time_steps;
    figure;
    plot(t, out.ms,'LineWidth',3);
    title('Ranks of approximations over time');
end
