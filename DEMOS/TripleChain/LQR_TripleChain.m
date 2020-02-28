function LQR_TripleChain(n1, usfs, shifts, istest)
% LQR_TripleChain computes the optimal feedback for the LQR problem (e.g.
% [1]) with respect to the Triple-Chain-Oscillator from [2]. The
% computations use either first or second companion form linearization. In
% order to reduce the complexity, the 'so_1' and 'so_2' function handles
% cast all operations in the linearized 2n x 2n model back to operations
% with respect ot the original nxn matrices [3,4]. The default function
% handles explicitly form the 2n x 2n matrices and treat them as a first
% order system.
%
% Usage:      LQR_TripleChain(n1, oper, shifts, istest)
% 
% Inputs: 
% 
% n1          length of a single chain in the model
%             (optional; defaults to 1000)
%
% usfs        the set of user supplied functions to use.
%             possible values: 'so_1', 'so_2', 'default'
%             (optional; defaults to 'so_1'
%
% shifts      the desired ADI shift selection strategy
%             possible values: 'heur', 'projection'
%             (optional; defaults to 'projection'
%
% istest      flag to determine whether this demo runs as a CI test or 
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% References:
%
% [1] A. Locatelli, Optimal Control: An Introduction, Birkhäuser, Basel,
%     Switzerland, 2001.
%
% [2] N. Truhar, K. Veselić, An efficient method for estimating the optimal
%     dampers’ viscosity for linear vibrating systems using Lyapunov
%     equation, SIAM J. Matrix Anal. Appl. 31 (1) (2009) 18–39.
%     https://doi.org/10.1137/070683052. 
%
% [3] P. Benner, J. Saak, Efficient Balancing based MOR for Second Order
%     Systems Arising in Control of Machine Tools, in: I. Troch, F.
%     Breitenecker (Eds.), Proceedings of the MathMod 2009, no. 35 in
%     ARGESIM-Reports, Vienna Univ. of Technology, ARGE Simulation News,
%     Vienna, Austria, 2009, pp. 1232–1243, iSBN/ISSN:978-3-901608-35-3. 
%
% [4] P. Benner, P. Kürschner, J. Saak, An improved numerical method for
%     balanced truncation for symmetric second order systems, Math. Comput.
%     Model. Dyn. Syst. 19 (6) (2013) 593–615.
%     https://doi.org/10.1080/13873954.2013.794363.   
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
narginchk(0,4);
if nargin<1, n1=1000; end
if nargin<2, usfs='so_1'; end
if nargin<3, shifts='projection'; end
if nargin<4, istest=0; end
%%
% set operation
oper = operatormanager(usfs);
% Problem data

alpha=2;
beta=5;
v=5;

switch usfs
    case 'so_1'
        [eqn.M_,eqn.E_,eqn.K_]=triplechain_MSD(n1,alpha,beta,v);
        s  = size(eqn.K_,1);
        eqn.B = [zeros(s,1); ones(size(eqn.K_,1),1)];
        eqn.C = [ones(1,size(eqn.K_,1)) zeros(1, s)];
    case 'so_2'
        [eqn.M_,eqn.E_,eqn.K_]=triplechain_MSD(n1,alpha,beta,v);
        s  = size(eqn.K_,1);
        eqn.B = [ ones(size(eqn.K_,1),1); zeros(s,1)];
        eqn.C = eqn.B';
    case 'default'
        [M_,E_,K_]=triplechain_MSD(n1,alpha,beta,v);
        s = size(K_,1);
        eqn.A_ = [sparse(s,s),-K_;-K_,-E_];
        eqn.E_ = [-K_,sparse(s,s);sparse(s,s),M_];
        eqn.B = [zeros(s,1); ones(size(K_,1),1)];
        eqn.C = [ones(1,size(K_,1)) zeros(1, s)];
        clear M_ E_ K_ s;
end

eqn.haveE=1;


%%
% ADI tolerances and maximum iteration number
opts.adi.maxiter = 200;
opts.adi.res_tol = 1e-10;
opts.adi.rel_diff_tol = 0;
opts.adi.info = 0;
opts.adi.accumulateK = 1;
opts.adi.accumulateDeltaK = 0;
opts.adi.compute_sol_fac = 1;

eqn.type = 'T';
%%
%Heuristic shift parameters via basic Arnoldi 
n=oper.size(eqn, opts);
switch shifts
    case 'heur'
        opts.shifts.num_desired=25;
        opts.shifts.num_Ritz=50;
        opts.shifts.num_hRitz=25;

        opts.shifts.info=0;
        opts.shifts.method = 'heur';
        opts.shifts.b0=ones(n,1);
    case 'projection'
        opts.shifts.num_desired=6;
        opts.shifts.method = 'projection';
        n=oper.size(eqn, opts);
        opts.shifts.b0=ones(n,1);
end
opts.shifts.truncate = 1e6; % remove all shifts larger than 1e6 or smaller 
                            % than 1e-6 in absolute value in order to avoid
                            % loosing information about M or K in the
                            % shifted coefficients (p^2*M-pD+K)  
%%
% Newton tolerances and maximum iteration number
opts.nm.maxiter = 25;
opts.nm.res_tol = 1e-10;
opts.nm.rel_diff_tol = 1e-16;
opts.nm.info = 1;
opts.nm.accumulateRes = 1;
opts.nm.linesearch = 1;
opts.nm.inexact = 'quadratic';
opts.nm.tau = 0.1;
opts.norm = 'fro';
%%
tic;
outnm = mess_lrnm(eqn, opts, oper);
toc;

if istest
    if min(outnm.res)>=opts.nm.res_tol
       error('MESS:TEST:accuracy','unexpectedly inaccurate result'); 
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


%% Lets try the RADI method and compare
% RADI-MESS settings
opts.shifts.history = opts.shifts.num_desired*size(eqn.C,1);

opts.shifts.method = 'projection';


% .. Suggest false (smart update is faster; convergence is the same).
opts.shifts.naive_update_mode = false; 
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
    if min(outradi.res)>=opts.radi.res_tol
       error('MESS:TEST:accuracy','unexpectedly inaccurate result'); 
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
if istest
    nrm = norm(outnm.K-outradi.K,'fro');
    nrmNM=norm(outnm.K,'fro');
    if nrm/nrmNM >= 1e-9
        error('MESS:TEST:accuracy',...
            'unexpectedly inaccurate result: ||K_NM - K_RADI||_F / ||K_NM||_F=%g',nrm/nrmNM);
    end 
else
    figure(3);
    ls_nm=[outnm.adi.niter];
    ls_radi=1:outradi.niter;

    semilogy(cumsum(ls_nm),outnm.res,'k--',ls_radi,outradi.res,'b-');
    title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
    xlabel('number of solves with A+p*M');
    ylabel('normalized residual norm');
    legend('LR-NM','RADI');
end