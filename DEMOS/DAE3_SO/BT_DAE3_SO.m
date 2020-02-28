function BT_DAE3_SO(model, tol, max_ord, maxiter, istest)
%%
% BT_DAE3_SO demonstrates the use of balanced truncation for index-3
% second order models. It reports the results presented in [3].
%
%   Usage:
%    BT_DAE3_SO(model,tol,max_ord,maxiter,test)
%
% Input:
%
%  model    choice of the example system. Possible values
%           'Stykel_small'   Stykel's mass spring damper system from [1]
%                            of dimension 600 in second order form (default) 
%           'Stykel_large'   Stykel's mass spring damper system from [1]
%                            of dimension 6000 in second order form 
%           'Truhar_Veselic' a version of the triple chain mass spring damper 
%                            test example from [2] with holonomic constraints
%                            coupling certain masses. Dimension 451 in second
%                            order form.
%           'TV'             another version of the example from [2] with
%                            other constraints
%           'TV2'            another version of the example from [2] with
%                            yet another set of constraints
%
%  tol       Balanced Truncation tolerance (optional, default: 1e-4)
%
%  max_ord   maximal reduced order allowed (optional, default: 500)
%
%  maxiter   maximal number of ADI iterations (optional, default: 300)
%
%  istest    decides whether the function runs as an interactive demo or a
%            continuous integration test. 
%            (optional; defaults to 0, i.e. interactive demo)
%
% References:
% [1] V. Mehrmann and T. Stykel. Balanced truncation model reduction for 
%     large-scale systems in descriptor form. 
%     In P. Benner, V. Mehrmann, and D. Sorensen, editors, Dimension
%     Reduction of Large-Scale Systems, volume 45 of Lecture Notes in 
%     Computational Science and Engineering, pages 83–115. Springer-Verlag, 
%     Berlin/Heidelberg, 2005.
%
% [2] N. Truhar and K. Veselic, An efficient method for estimating the 
%     optimal dampers’ viscosity for linear vibrating systems using 
%     Lyapunov equation, SIAM J. Matrix Anal. Appl., 31 (2009), pp. 18–39.
%
% [3] J. Saak, M. Voigt, Model reduction of constrained mechanical systems
%     in M-M.E.S.S., IFAC-PapersOnLine 9th Vienna International Conference
%     on Mathematical Modelling MATHMOD 2018, Vienna, Austria, 21–23
%     February 2018 51 (2) (2018) 661–666.
%     https://doi.org/10.1016/j.ifacol.2018.03.112.      

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

%% Input checks
narginchk(0,5);
if nargin==0
    model='Stykel_small';
end

if nargin < 2
    tol = 1e-4; 
end

if nargin < 3
    max_ord = 500;
end

if nargin < 4
    maxiter = 300;
end

if nargin < 5
    istest = 0; 
end
%% set operation manager for the structured computations of Gramians
oper = operatormanager('dae_3_so');

%% load problem data

switch lower(model)
    case {'stykel_small','stykel_large'}
        if strcmp(model,'stykel_small')
            sys = load(sprintf('%s/../models/ms_ind3_by_t_stykel/g600.mat',...
                fileparts(mfilename('fullpath'))));
        else
            sys = load(sprintf('%s/../models/ms_ind3_by_t_stykel/g6000.mat',...
                fileparts(mfilename('fullpath'))));
        end
        eqn.M_=sys.M;
        eqn.E_=sys.D;
        eqn.K_=sys.K;
        eqn.G_=sys.G;
        eqn.haveE=1;
        eqn.alpha = -0.02;
        nv = size(eqn.M_,1);
        np = size(eqn.G_,1);
        eqn.B = full(sys.B(1:2*nv,:));
        eqn.C = full(sys.C(:,1:2*nv));
        clear E A B C M D K G;        
    case 'tv2'
        n1=151; % make sure it is odd!
        alpha=0.01;
        beta=alpha;
        v=5e0;
        
        [eqn.M_,eqn.E_,eqn.K_]=triplechain_MSD(n1,alpha,beta,v);
        eqn.E_ = -eqn.E_;
        eqn.K_ = -eqn.K_;
        eqn.G_ = sparse(3,3*n1+1); 
        n12=ceil(n1/2); % n1 is odd so this is the index of 
                        % the center of the string
        eqn.G_(1,1)=-1;       eqn.G_(1,n12)=2;       eqn.G_(1,n1)=-1;
        eqn.G_(2,n1+1)=-1;    eqn.G_(2,n1+n12)=-1;   eqn.G_(2,2*n1)=-1;
        eqn.G_(3,2*n1+1)=-1;  eqn.G_(3,2*n1+n12)=2;  eqn.G_(3,3*n1)=-1;
       
        nv = size(eqn.M_,1);
        np = size(eqn.G_,1);
        eqn.B=zeros(6*n1+2,3);
        eqn.B(3*n1+6,1)=1;
        eqn.B(5*n1+n12+6,2)=1;
        eqn.B(end-6,3)=1;
        eqn.C=zeros(3,6*n1+2);
        eqn.C(1,3*n1+1)=1;
        eqn.C(2,2*n1)=1;
        eqn.C(3,2*n1+floor(n12/2))=1;
        eqn.haveE=1;
        eqn.alpha = -0.02;
    case 'tv'
        n1=151; % make sure it is odd!
        alpha=0.01;
        beta=alpha;
        v=5e0;
        
        [eqn.M_,eqn.E_,eqn.K_]=triplechain_MSD(n1,alpha,beta,v);
        eqn.E_ = -eqn.E_;
        eqn.K_ = -eqn.K_;
        eqn.G_ = sparse(3,3*n1+1); 
        n12=ceil(n1/2); % n1 is odd so this is the index of
                        % the center of the string
        eqn.G_(1,1)=1/3;       eqn.G_(1,n12)=1/3;       eqn.G_(1,n1)=1/3;
        eqn.G_(2,n1+1)=1/3;    eqn.G_(2,n1+n12)=1/3;   eqn.G_(2,2*n1)=1/3;
        eqn.G_(3,2*n1+1)=1/3;  eqn.G_(3,2*n1+n12)=1/3;  eqn.G_(3,3*n1)=1/3;
       
        nv = size(eqn.M_,1);
        np = size(eqn.G_,1);
        eqn.B=[zeros(3*n1+1,1);ones(3*n1+1,1)];
        eqn.C=[ones(3*n1+1,1);zeros(3*n1+1,1)]';
        eqn.haveE=1;
        eqn.alpha = -0.02;
   case 'truhar_veselic'
        n1=1500; % make sure it is even!
        alpha=0.01;
        beta=alpha;
        v=5e0;
             
        [eqn.M_,eqn.E_,eqn.K_]=triplechain_MSD(n1,alpha,beta,v);
        eqn.E_ = -eqn.E_;
        eqn.K_ = -eqn.K_;
        eqn.G_ = sparse(6,3*n1+1); 
        eqn.G_(1,1)=1;   eqn.G_(1,n1/2)=-1;
        eqn.G_(2,n1/2+1)=1;eqn.G_(2,n1)=-1;
        eqn.G_(3,n1+1)=1;eqn.G_(3,n1+n1/2)=-1;
        eqn.G_(4,n1+n1/2+1)=1;eqn.G_(4,2*n1)=-1;
        eqn.G_(5,2*n1+1)=1;eqn.G_(5,2*n1+n1/2)=-1;
        eqn.G_(6,2*n1+n1/2+1)=1;eqn.G_(6,3*n1)=-1;
        nv = size(eqn.M_,1);
        np = size(eqn.G_,1);
        eqn.B=[zeros(3*n1+1,1);ones(3*n1+1,1)];
        eqn.C=[ones(3*n1+1,1);zeros(3*n1+1,1)]';
        eqn.haveE=1;
        eqn.alpha = -0.02;
    otherwise
        fprintf('unknown model requested!\n');
        return
end
%% options
% Adi options
opts.adi.maxiter = maxiter;
opts.adi.res_tol = 1e-10;
opts.adi.rel_diff_tol = 1e-11;
opts.norm = 'fro';
opts.shifts.method='projection';
opts.shifts.num_desired=25;

%% LRADI for the two Gramian factors
%  controllability Gramian
eqn.type = 'T';
opts.adi.info = 1;
tic;
opts.shifts.p = mess_para(eqn, opts, oper);
% use an additional alpha-shift to improve convergence and ROM quality for
% the triple chain model
if strcmp(model,'Truhar_Veselic')||strcmp(model,'TV')||strcmp(model,'TV2')
    opts.shifts.p= opts.shifts.p-0.5;
end

outC = mess_lradi(eqn, opts, oper);
toc;

% observability Gramian
eqn.type = 'N';
tic;
outB = mess_lradi(eqn, opts, oper);
toc;

%% Reduced Order Model computation via square root method (SRM)
fprintf('\nComputing reduced order model via square root method\n\n');
opts.srm.tol=tol;
opts.srm.max_ord=max_ord;
opts.srm.info=1;

tic;
[TL, TR] = mess_square_root_method(eqn, opts, oper, outB.Z, outC.Z);
% compute ROM matrices
ROM.A = TL' * oper.mul_A(eqn, opts, 'N', TR, 'N');
ROM.B = TL' * eqn.B;
ROM.C = eqn.C * TR;
ROM.E = eye(size(ROM.A));
toc;

%% Frequency-domain evaluation of the (transfer function of the)
%  ROM and comparison to the original model. 
%
% We feed the mess_sigma_plot with usfs that do not exploit the DAE structure: 
tic;
opts.sigma.nsample = 200;
if istest
    opts.sigma.info = 0;
else
    opts.sigma.info = 2;
end
if strcmp(model,'TV2')
    opts.sigma.fmin=-2;
    opts.sigma.fmax=1;
else
    opts.sigma.fmin=-4;
    opts.sigma.fmax=4;
end

NG = sparse(np,nv);
NS = sparse(np,np);
eqnu.M_ = [eqn.M_ NG'; NG NS];
eqnu.E_ = [-eqn.E_ NG'; NG NS];
eqnu.K_ = [-eqn.K_ -eqn.G_';-eqn.G_ NS];
eqnu.C =  [zeros(size(eqn.C,1),np+nv), eqn.C(:,1:nv), zeros(size(eqn.C,1),np)];
eqnu.B =  [eqn.B(nv+1:2*nv,:);zeros(2*np+nv,size(eqn.B,2))];
eqnu.haveE = 1;

operu = operatormanager('so_1');

out = mess_sigma_plot(eqnu, opts, operu, ROM); err = out.err;

toc;
%% final accuracy test used in the continuous integration system or
%  plot of the computed  
if istest
    % the errors are not always perfect in this example, but let's see
    % wether they are "good enough"...
    if (max(err) > 50*tol)
        error('MESS:TEST:accuracy', ['unexpectedly inaccurate result ' ...
                            'for %s %g %d %d (%g)'], model, tol, ...
              max_ord, maxiter,max(err));  
    end
end