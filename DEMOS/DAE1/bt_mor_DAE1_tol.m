function bt_mor_DAE1_tol(istest)
% Computes a reduced order model via Balanced Truncation for the proper
% index-1 System BIPS98_606 from https://sites.google.com/site/rommes/software
% following the method suggested in [1]
%
% Input: 
% istest    decides whether the function runs as an interactive demo or a
%           continuous integration test. (optional; defaults to 0, i.e.
%           interactive demo)
%
% References:
%[1] F. Freitas, J. Rommes, N. Martins, Gramian-based reduction method
%    applied to large sparse power system descriptor models, IEEE Trans.
%    Power Syst. 23 (3) (2008) 1258–1270. doi:10.1109/TPWRS.2008.926693. 

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

%%
if nargin<1, istest=0; end
%%
% set operation manager for the Gramian computations
oper = operatormanager('dae_1');

%% Read problem data
fname  = sprintf('%s/../models/BIPS/bips98_606.mat',...
    fileparts(mfilename('fullpath')));
Bips = load(fname);
% from https://sites.google.com/site/rommes/software
% Note that we here use the alpha shift suggested by Rommes and coauthors
% [1], but do not perform the shift back.
p = find(diag(Bips.E));
np = find(diag(Bips.E) == 0);
pp = [p;np];
eqn.A_ = Bips.A(pp, pp)-0.05*Bips.E(pp,pp);
eqn.E_ = Bips.E(pp, pp);
eqn.B = Bips.b(pp, :);
eqn.C = Bips.c( : , pp);
eqn.st = length(p);
eqn.haveE = 1;
%% Turn off  close to singular warnings
% (this model is really badly conditioned)
orig_warnstate = warning('OFF','MATLAB:nearlySingularMatrix');

%%
% ADI tolerances and maximum iteration number
opts.adi.maxiter = 300;
opts.adi.res_tol = 1e-10;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 1;
opts.norm = 'fro';

%%
opts.shifts.method = 'projection';
opts.shifts.num_desired= 20;

opts.shifts.p=mess_para(eqn,opts,oper);

disp(opts.shifts.p);

%%
eqn.type='N';
tic;
outB = mess_lradi(eqn, opts, oper);
toc;

if istest
    if min(outB.res)>=opts.adi.res_tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result'); 
    end
else
    figure;
    semilogy(outB.res);
    title('0= BB^T + AXM^T + MXA^T');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
end

disp('size outB.Z:');
disp(size(outB.Z));

%%
eqn.type='T';
tic;
outC = mess_lradi(eqn, opts, oper);
toc;


if istest
    if min(outC.res)>=opts.adi.res_tol
        error('MESS:TEST:accuracy','unexpectedly innacurate result'); 
    end
else
    figure;
    semilogy(outC.res);
    title('0= C^TC + A^TXE + E^TXA');
    xlabel('number of iterations');
    ylabel('normalized residual norm');
    pause(1);
end
disp('size outC.Z:');
disp(size(outC.Z));

%% Compute reduced system matrices
% Perform Square Root Method  (SRM)

% BT tolerance and maximum order for the ROM
opts.srm.tol=1e-3;
opts.srm.max_ord=250;

% SRM verbosity
if istest
    opts.srm.info=1;
else
    opts.srm.info=2;
end

%The actual SRM
[TL,TR,hsv] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);

% compute ROM matrices
B1 = TL'*(eqn.A_(1:eqn.st,1:eqn.st))*TR;
B2 = TL'*(eqn.A_(1:eqn.st,eqn.st+1:end));
A1 = eqn.A_(eqn.st+1:end,1:eqn.st)*TR;

ROM.A = B1 - B2*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\A1);
ROM.B = TL'*eqn.B(1:eqn.st,:) - ...
     B2*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\eqn.B(eqn.st+1:end,:));
ROM.C = eqn.C(:,1:eqn.st)*TR - ...
     eqn.C(:,eqn.st+1:end)*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\A1);
ROM.D = -eqn.C(:,eqn.st+1:end)*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\...
    eqn.B(eqn.st+1:end,:));
ROM.E = eye(size(ROM.A));
%% Evaluate the ROM quality
% while the Gramians are computed on the hidden manifold, we need to do the
% frequency domain computations without (implicitly) using the Schur
% complement (due to the construction of the function handles)
oper = operatormanager('default');

if istest
    opts.sigma.info=0;
else
    opts.sigma.info=2;
end

opts.sigma.fmin=-3;
opts.sigma.fmax=4;

err = mess_sigma_plot(eqn, opts, oper, ROM);

if istest
    if max(err)>5e-3
        error('MESS:TEST:accuracy','unexpectedly innacurate result'); 
    end   
else
    figure;
    semilogy(hsv);
    title('Computed Hankel singular values');
    xlabel('index');
    ylabel('magnitude');
end
%% reset warning state
warning(orig_warnstate);