function out =  LQR_rail_splitting(k, exp_action, method)

% Computes the optimal feedback via low-rank splitting schemes [1, 2] for 
% the selective cooling of Steel profiles application described in [3,4,5].
% 
% Inputs: 
% 
% k           refinement level of the model to use 
%             (1-4, i.e. 1357-79841 DOFs)
%             (optional, defaults to 1)
%
% exp_action  parameters for computing actions of matrix exponentials, see
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
%     (2015), pp. 2791-2796. https://doi.org/10.1109/TAC.2015.2398889.
% 
% [2] T. Stillfjord, Adaptive high-order splitting schemes for large-scale
%     diffederential Riccati equations, Numer. Algorithms,  (2017).
%     https://doi.org/10.1007/s11075-017-0416-8.
%
% [3] J. Saak, Effiziente numerische L\"{o}sung eines
%     Optimalsteuerungsproblems f\"{u}r die Abk\"{u}hlung von 
%     Stahlprofilen, Diplomarbeit, Fachbereich 3/Mathematik und Informatik, 
%     Universit\"{a}t Bremen, D-28334 Bremen (Sep. 2003).   
%
% [4] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353-356. https://doi.org/10.1007/3-540-27909-1_19. 
%
% [5] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universit\"{a}t Chemnitz, Chemnitz, Germany (Jul. 2009).  
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

if nargin < 1
    k = 1;
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

%% Equation parameters
% Default (E, A, B, C) system
oper = operatormanager('default');

eqn = getrail(k);
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
tic;
[out, ~, opts, ~] = mess_splitting_dre(eqn, opts, oper);
toc;

%%
t = opts.splitting.time_steps;
figure;
plot(t, out.ms);
title('Ranks of approximations over time');