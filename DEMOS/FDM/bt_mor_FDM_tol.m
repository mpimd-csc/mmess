function [Ar, Br, Cr] = bt_mor_FDM_tol(tol,n0,shifts,istest)
% bt_mor_FDM_tol computes a reduced order model via the standard Lyapunov
% balanced truncation (see e.g. [1]) for a finite difference discretized  
% convection diffusion model on the unit square described in [2].
%
% Usage: 
%    [Ar, Br, Cr] = bt_mor_FDM_tol(tol,max_ord,n0,test)
%
% Inputs
%
% tol         truncation tolerance for the Hankel singular values
%             (optional; defalts to 1e-6)
% 
% n0          n0^2 gives the dimension of the original model, i.e. n0 is
%             the number of degrees of freedom, i.e. grid points, per
%             spatial direction 
%             (optional; defaults to 50)
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
% [2] T. Penzl, Lyapack Users Guide, Tech. Rep. SFB393/00-33, 
%     Sonderforschungsbereich 393 Numerische Simulation auf massiv 
%     parallelen Rechnern, TU Chemnitz, 09107 Chemnitz, Germany, 
%     available from http://www.tu-chemnitz.de/sfb393/sfb00pr.html. (2000).


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
% BT tolerance and maximum order for the ROM
if nargin<1, tol=1e-6; end
if nargin<2, n0 = 50; end
if nargin<3, shifts='heur'; end
if nargin<4, istest = 0; end

% ADI tolerance and maximum iteration number
opts.adi.maxiter = 100;
opts.adi.res_tol = 1e-9;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 1;
opts.norm = 'fro';

%%
% operations
oper = operatormanager('default');

% Problem data
eqn.A_ = fdm_2d_matrix(n0,'10*x','100*y','0');
eqn.B = fdm_2d_vector(n0,'.1<x<=.3');
eqn.C = fdm_2d_vector(n0,'.7<x<=.9');eqn.C=eqn.C';
eqn.haveE = 0;

% problem dimension
n=oper.size(eqn, opts);

%%
%Heuristic Parameters via basic Arnoldi or Projection shifts
opts.shifts.method = shifts;
switch opts.shifts.method
    case 'heur'
        opts.shifts.num_desired=25;
        opts.shifts.num_Ritz=50;
        opts.shifts.num_hRitz=25;
        opts.shifts.b0=ones(n,1);
    case 'projection'
        opts.shifts.method = 'projection';
        opts.shifts.num_desired = 10;
    otherwise
        warning('MESS:Test:invalid_parameter_fallback',...
            'invalid shift selection, falling back to "heur"');
        opts.shifts.num_desired=25;
        opts.shifts.num_Ritz=50;
        opts.shifts.num_hRitz=25;
        opts.shifts.b0=ones(n,1);
        opts.shifts.method = 'heur';
end

opts.shifts.p=mess_para(eqn,opts,oper);

disp('(Initial) ADI shifts');
disp(opts.shifts.p);
%%
%observability
eqn.type='N';
tic;
outB  = mess_lradi(eqn, opts, oper);
toc;

if istest
    if min(outB.res)>=opts.adi.res_tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result');
    end
else
    figure(1);
    semilogy(outB.res);
    title('AX + XA^T = -BB^T');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end

disp('size outB.Z:');
disp(size(outB.Z));

%%
%controllability
eqn.type = 'T';
tic;
outC = mess_lradi(eqn, opts, oper);
toc;

if istest
    if min(outC.res)>=opts.adi.res_tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result');
    end
else
    figure(2);
    semilogy(outC.res);
    title('A^TX + XA = -C^TC');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outC.Z:');
disp(size(outC.Z));

%%
opts.srm.tol=tol;
opts.srm.info=1;
[TL, TR, hsv, eqn, opts, ~] = mess_square_root_method(eqn,opts,oper,...
    outB.Z,outC.Z);

Ar = TL'*(eqn.A_*TR);
Br = TL'*eqn.B;
Cr = eqn.C*TR;

opts.sigma.nsample = 200;  % 200 frequency samples
opts.sigma.fmin = -3;      % min. frequency 1e-3
opts.sigma.fmax = 4;       % max. frequency 1e4
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