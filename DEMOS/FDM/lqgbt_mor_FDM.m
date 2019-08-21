function [Ar, Br, Cr] = lqgbt_mor_FDM(tol,max_ord,n0,istest)
% LQGBT_MOR_FDM computes a reduced order model via the lienar quadratic 
% Gaussian balanced truncation [1] for a finite difference discretized 
% convection diffusion model on the unit square described in [2].
%
% Usage: 
%    [Ar, Br, Cr] = lqgbt_mor_FDM(tol,max_ord,n0,test)
%
% Inputs
%
% tol         truncation tolerance for the LQG characteristic values
% 
% max_ord     maximum allowed order for the reuced order model
%
% n0          n0^2 gives the dimension of the original model, i.e. n0 is
%             the number of degrees of freedom per spatial direction
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
% [1] D. Mustafa, K. Glover, Controller design by H∞ -balanced truncation,
%     IEEE Trans. Autom. Control 36 (6) (1991) 668–682. 
%     https://doi.org/10.1109/9.86941. 
%
% [2] T. Penzl, Lyapack Users Guide, Tech. Rep. SFB393/00-33, 
%     Sonderforschungsbereich 393 Numerische Simulation auf massiv 
%     parallelen Rechnern, TU Chemnitz, 09107 Chemnitz, Germany, 
%     available from http://www.tu-chemnitz.de/sfb393/sfb00pr.html. (2000).

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
%               2009-2017
%

%%
narginchk(0,4);
%%
% LQGBT Tolerance and maximum order of the ROM.
if nargin<1
    tol     = 1e-6;
end
if nargin<2
    max_ord = 250;
end
if nargin<3
    n0 = 60; % n0 = number of grid points in either space direction; 
             % n = n0^2 is the problem dimension!
             % (Change n0 to generate problems of different size.)
end
if nargin<4
    istest=0;
end
% Problem data

eqn.A_ = fdm_2d_matrix(n0,'10*x','100*y','0');
eqn.B  = fdm_2d_vector(n0,'.1<x<=.3');
eqn.C  = [fdm_2d_vector(n0,'.7<x<=.8'), fdm_2d_vector(n0,'.8<x<=.9')];
eqn.C  = eqn.C';

%%
% Set operations.
oper = operatormanager('default');

%%
% ADI options.
opts.adi.maxiter          = 200;
opts.adi.res_tol           = 1.0e-10;
opts.adi.rel_diff_tol            = 0;
opts.adi.info             = 0;
opts.adi.accumulateK      = 1;
opts.adi.accumulateDeltaK = 0;
opts.adi.compute_sol_fac  = 1;

n = oper.size(eqn, opts);
opts.shifts.num_desired = 25;
opts.shifts.num_Ritz = 50;
opts.shifts.num_hRitz = 25;

opts.shifts.info   = 0;
opts.shifts.method = 'heur';
opts.shifts.b0     = ones(n,1);

%%
opts.norm          = 'fro';
% Newton options.
opts.nm.maxiter       = 25;
opts.nm.res_tol        = 1e-10;
opts.nm.rel_diff_tol         = 1e-16;
opts.nm.info          = 1;
opts.nm.accumulateRes = 1;
opts.nm.linesearch    = 1;
opts.nm.inexact       = 'quadratic';
opts.nm.tau           = 0.1;

%%
% Solve the filter Riccati equation.
%   A*X + X*A' - X*C'*C*X + B*B' = 0
tic;
eqn.type = 'N';
outC = mess_lrnm(eqn, opts, oper);
toc;

if istest
    if min(outC.res)>=opts.nm.res_tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result');
    end
else
    figure(1);
    semilogy(outC.res);
    title('AX + XA^T - XC^TCX + BB^T = 0');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end


%%
% Solve the regulator Riccati equation.
%   A'*X + X*A - X*B*B'*X + C'*C = 0
tic;
eqn.type = 'T';
outB = mess_lrnm(eqn, opts, oper);
toc;

if istest
    if min(outB.res)>=opts.nm.res_tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result');
    end
else
    figure(2);
    semilogy(outB.res);
    title('A^TX + XA - XBB^TX + C^TC = 0');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end


%%
% % Model reduction by square root method.
opts.srm.tol=tol;
opts.srm.max_ord = max_ord; 
opts.srm.info=1;
[TL, TR, hsv, eqn, opts, ~] = mess_square_root_method(eqn,opts,oper,...
    outB.Z,outC.Z);

Ar = TR'*(eqn.A_*TL);
Br = TR'*eqn.B;
Cr = eqn.C*TL;

opts.sigma.nsample = 200;  % 200 frequency samples
opts.sigma.fmin = -2;      % min. frequency 1e-3
opts.sigma.fmax = 6;       % max. frequency 1e4
if istest
    opts.sigma.info=1;     % no output
else
    opts.sigma.info = 2;   % show messages and plots
end
ROM.A = Ar;
ROM.B = Br;
ROM.C = Cr;
ROM.E = eye(size(ROM.A,1));

err = mess_sigma_plot(eqn, opts, oper, ROM);

%%
% Report.


if istest
    if max(err)>=tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result');
    end    
end
