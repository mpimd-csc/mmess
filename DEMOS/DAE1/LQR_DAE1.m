function LQR_DAE1(istest)
% Computes a Riccati feedback control for the proper
% index-1 System BIPS98_606 from https://sites.google.com/site/rommes/software
% following the ideas introduced in [1] for Lyapunov equations using he
% Newton-ADI iteration.
%
% Input: 
% istest  decides whether the function runs as an interactive demo or a
%         continuous integration test. (optional; defaults to 0, i.e.
%         interactive demo)
%
% References:
%[1] F. Freitas, J. Rommes, N. Martins, Gramian-based reduction method
%    applied to large sparse power system descriptor models, IEEE Trans.
%    Power Syst. 23 (3) (2008) 1258–1270. doi:10.1109/TPWRS.2008.926693. 

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
if nargin<1, istest=0; end

%%
% set operation
oper = operatormanager('dae_1');

%% Problem data
fname  = sprintf('%s/../models/BIPS/bips98_606.mat',...
    fileparts(mfilename('fullpath')));
Bips = load(fname);
% from https://sites.google.com/site/rommes/software
p = find(diag(Bips.E));
np = find(diag(Bips.E) == 0);
pp = [p;np];
eqn.A_ = Bips.A(pp, pp);
eqn.E_ = Bips.E(pp, pp);
eqn.B = Bips.b(pp, :);
eqn.C = 0.01*Bips.c( : , pp);
eqn.st = length(p);
eqn.haveE = 1;
clear Bips;
%% Turn off  close to singular warnings
%  (this model is really badly conditioned)
orig_warnstate = warning('OFF','MATLAB:nearlySingularMatrix');

%%
opts.norm = 'fro';

% ADI tolerances and maximum iteration number
opts.adi.maxiter = 300;
opts.adi.res_tol = 1e-12;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;
opts.adi.projection.freq=0;

eqn.type='N';
%%
opts.shifts.num_desired=25;
opts.shifts.num_Ritz=50;
opts.shifts.num_hRitz=25;
opts.shifts.method = 'projection';
opts.shifts.num_desired= 9;

%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 30;
opts.nm.res_tol = 1e-10;
opts.nm.rel_diff_tol = 1e-16;
opts.nm.info = 1;
opts.nm.linesearch = 1;
opts.nm.accumulateRes = 1;
%%
tic;
outnm = mess_lrnm(eqn, opts, oper);
toc;
if istest
    if min(outnm.res)>=opts.nm.res_tol
        error('MESS:TEST:accuracy','unexpectedly inaccurate result in LRNM'); 
    end
else
    figure(1);
    semilogy(outnm.res);
    title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outnm.Z:');
disp(size(outnm.Z));

%% Lets try RADI
opts.norm = 2;

% RADI-MESS settings
opts.shifts.history = opts.shifts.num_desired*size(eqn.C,1);
opts.shifts.method  = 'gen-ham-opti';

opts.shifts.naive_update_mode = false; 
% .. Suggest false (smart update is faster; convergence is the same).

opts.radi.compute_sol_fac = 1;
opts.radi.get_ZZt = 1;
opts.radi.compute_res = 0;
opts.radi.maxiter = 500;
opts.radi.res_tol = opts.nm.res_tol;
opts.radi.rel_diff_tol = 0;
opts.radi.info = 1;

tic;
outradi = mess_lrradi( eqn, opts, oper );
toc;

if istest
    if min(outradi.res)>=opts.radi.res_tol
        error('MESS:TEST:accuracy','unexpectedly inaccurate result in RADI'); 
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

%% reset warning state
warning(orig_warnstate');
