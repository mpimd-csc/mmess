function bt_mor_DAE2(problem,lvl,re,istest)
% Computes a standard ROM by implicitly solving the generalized Lyapunov
% equations for the equivalent projected system on the hidden manifold.
%
% Inputs:
% problem       either 'Stokes' or 'NSE' to choose the Stokes demo or the
%               linearized Navier-Stokes-Equation. (required)
%
% lvl           discretization level 1 through 5 
%               (optional, only used in 'NSE' case, default: 1)
%
% re            Reynolds number 300, 400, or 500
%               (optional, only used in 'NSE' case, default: 500)
%
% istest        flag to determine whether this demo runs as a CI test or 
%               interactive demo
%               (optional, defaults to 0, i.e. interactive demo)
% 
% Note that the 'NSE' option requires additional data available in a
% separate 270MB archive and at least the 5th discretization level needs a
% considerable amount of main memory installed in your machine.
%
% See 
% P. Benner, J. Saak, M. M. Uddin, Balancing based model reduction for 
% structured index-2 unstable descriptor systems with application to flow
% control, Numerical Algebra, Control and Optimization 6 (1) (2016) 1–20. 
% https://doir.org/10.3934/naco.2016.6.1.


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

% ADI tolerance and maximum iteration number
opts.adi.maxiter = 350;
opts.adi.res_tol = sqrt(eps);
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 1;
opts.shifts.info = 1;
opts.norm = 'fro';

oper = operatormanager('dae_2');
%% Problem data
if nargin<1, problem='stokes'; end
if nargin<2, lvl=1; end
if nargin<3, re=500; end
if nargin<4, istest=0; end

switch lower(problem)
    case 'stokes'
        nin = 5;
        nout = 5;
        nx = 10;
        ny = 10;
        [eqn.E_,eqn.A_,eqn.Borig,eqn.Corig]=stokes_ind2(nin,nout,nx,ny);
        n=size(eqn.E_,1);
        eqn.haveE=1;
        st=full(sum(diag(eqn.E_)));
        eqn.st=st;
        eqn.B=eqn.Borig(1:st,:);
        eqn.C=eqn.Corig(:,1:st);
    case 'nse'
        try
            load(sprintf('%s/../models/NSE/mat_nse_re_%d',...
                    fileparts(mfilename('fullpath')),re),'mat');
        catch
            error(['The files mat_nse_re_300.mat, mat_nse_re_400.mat ', ...
                'and mat_nse_re_500.mat are available for dowload in ', ...
                'a separate archive (270MB each). Please fetch them ', ...
                'from the MESS download page and unpack them into ', ...
                'the DEMOS/models/NSE folder.']);
        end
        eqn.A_=mat.mat_v.fullA{lvl};
        eqn.E_=mat.mat_v.E{lvl};
        eqn.B=mat.mat_v.B{lvl};
        eqn.C=mat.mat_v.C{lvl};
        eqn.st=mat.mat_mg.nv(lvl);
        st=eqn.st;
        eqn.haveE=1;
        n=size(eqn.E_,1);
end
%%
eqn.type='N';
eqn.G=eqn.B;
if strcmp(problem,'NSE') && (re>200)
    eqn.V=-mat.mat_v.Feed_0{lvl};
    eqn.U = eqn.B;
    eqn.haveUV=1;
end


opts.shifts.num_desired=6;
opts.shifts.num_Ritz=40;
opts.shifts.num_hRitz=40;
opts.shifts.method='projection';

opts.shifts.b0=ones(size(eqn.A_,1),1);

opts.shifts.p=mess_para(eqn,opts,oper);

tic;
outB = mess_lradi(eqn,opts,oper);
toc;

if not(istest)
    figure(1);
    semilogy(outB.res);
    title('AXM^T + MXA^T = -BB^T');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outB.Z:');
disp(size(outB.Z));

%%
eqn.type = 'T';
eqn.G=eqn.C';
if strcmp(problem,'NSE') && (re>200)
    if (re == 500) 
        % the dataset stores the feedback of the adjoint system transposed
        % for Reynolds 500
        eqn.U=-mat.mat_v.Feed_1{lvl}';
    else
        eqn.U=-mat.mat_v.Feed_1{lvl};
    end
    eqn.V = eqn.C';
    eqn.haveUV=1;
end


opts.shifts.num_desired=6;
opts.shifts.num_Ritz=40;
opts.shifts.num_hRitz=40;
opts.shifts.method='projection';

opts.shifts.b0=ones(size(eqn.A_,1),1);

opts.shifts.p=mess_para(eqn,opts,oper);

tic;
outC = mess_lradi(eqn, opts, oper);
toc;

if not(istest)
    figure(2);
    semilogy(outC.res);
    title('A^TXM + M^TXA = -C^TC');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outC.Z:');
disp(size(outC.Z));

%% Compute reduced system matrices
% Perform Square Root Method  (SRM)

% BT tolerance and maximum order for the ROM
tic;
opts.srm.tol=1e-5;
opts.srm.max_ord=250;

% SRM verbosity
if istest
    opts.srm.info=1;
else
    opts.srm.info=2;
end

%The actual SRM
[TL,TR,hsv] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);

%%
ROM.A = TL'*(eqn.A_(1:st,1:st)*TR);
ROM.B = TL'*eqn.B(1:st,:);
ROM.C = eqn.C(:,1:st)*TR;

toc;
%%
tic;
%% Evaluate the ROM quality
% while the Gramians are computed exploiting the DAE structure, due to the
% construction of the function handles we can not do so for the transfer
% function. Therfore we need to extend the matrices B and C and call the
% 'default' usfs for unstructured computation:
switch lower(problem)
    case 'stokes'
        eqn.B=eqn.Borig;
        eqn.C=eqn.Corig;
    case 'nse'
        n = size(eqn.A_,1);
        eqn.B(st+1:n,:) = zeros(n-st,size(eqn.B,2));
        eqn.C(:,st+1:n) = zeros(size(eqn.C,1),n-st);
end
oper = operatormanager('default');

if istest
    opts.sigma.info=0;
else
    opts.sigma.info=2;
end

opts.sigma.fmin=-3;
opts.sigma.fmax=4;

err = mess_sigma_plot(eqn, opts, oper, ROM);

toc;
%%  
if istest
    if max(err)>=opts.srm.tol, error('MESS:TEST:accuracy','unexpectedly innacurate result'); end
else
    figure;
    semilogy(hsv);
    title('Computed Hankel singular values');
    xlabel('index');
    ylabel('magnitude');
end    
%%
fprintf(['\nComputing open loop step response of original and reduced order ' ...
    'systems and time domain MOR errors\n']);
open_step(eqn,ROM.A,ROM.B,ROM.C,problem,istest);
%%
fprintf('\nComputing ROM based feedback\n');
if exist('care', 'file')
    [~,~,Kr]=care(ROM.A,ROM.B,ROM.C'*ROM.C,eye(size(ROM.B,2)));
else
    Y = care_nwt_fac([],ROM.A,ROM.B,ROM.C,1e-12,50);
    Kr = (Y*ROM.B)'*Y;
end
K=[Kr*TL'*eqn.E_(1:st,1:st),zeros(size(Kr,1),n-st)];
%%
fprintf(['\nComputing closed loop step response of original and reduced order ' ...
    'systems and time domain MOR errors\n']);
closed_step(eqn,ROM.A,ROM.B,ROM.C,problem,K,Kr,istest);

