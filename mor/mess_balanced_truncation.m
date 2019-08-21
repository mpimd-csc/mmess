function [Er,Ar,Br,Cr,outinfo]=mess_balanced_truncation(E,A,B,C,max_order,trunc_tol,info,opts)
% Lyapunov Balanced truncation for descriptor systems with invertible E.
%
%  [Er,Ar,Br,Cr,outinfo]=mess_balanced_truncation(E,A,B,C,max_order,trunc_tol,info,opts)
%
% INPUTS:
%  E,A,B,C    The mass, system, input and output matrices describing the
%             original system
%  max_ord    maximum reduced order allowed 
%             (optional, defaults to size(A,1))
%  trunc_tol  error tolerance used for the Hankel singular value truncation
%             (optional, defaults to 1e-5)
%  info       verbosity control parameter (optional):
%             0  quiet (default)
%             1  show iteration numbers and residuals
%             >1 plot residual history
%             >2 compute and show the sigma and error plots
%  opts       options structure that can be used to pass setting to the
%             LRADI, ADI shift computation, or the square root method (optional)
%             (see corresponding routines for additional information)
%             It also has fields to control the plotting in case info>2:
%             opts.sigma.fmin      minimum value in the logspace for the 
%                                  sigma and error plots
%             opts.sigma.fmax      maximum value in the logspace for the 
%                                  sigma and error plots
%             opts.sigma.nsample   number of elements in the losgspace.
%
% OUTPUTS:
% Er, Ar, Br, Cr         the reduced order model matrices
% outinfo.TL outinfo.TR  the left and right transformation matrices
% outinfo.errbound       H-infinity error bound
% outinfo.hsv            vector with the computed Hankel singular values
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

narginchk(4,8);
%%
% control verbosity of the computations
if nargin<7
  opts.adi.info=0;
  opts.bt.info=0;
  opts.srm.info=0;
  info=0;
  opts.shifts.info = 0;
else
    if nargin <8
        opts.adi.info=info;
        opts.bt.info=info;
        opts.srm.info=info;
        opts.shifts.info=info;
    else
        if not(isfield(opts,'adi')) || not(isfield(opts.adi,'info'))
            opts.adi.info = info;
        end
        if not(isfield(opts,'bt')) || not(isfield(opts.bt,'info'))
            opts.bt.info = info;
        end
        if not(isfield(opts,'srm')) || not(isfield(opts.srm,'info'))
            opts.srm.info = info;
        end
        if not(isfield(opts,'shifts')) || not(isfield(opts.shifts,'info'))
            opts.shifts.info = info;
        end
    end
end

% BT tolerance and maximum order for the ROM
if (nargin<6 || isempty(trunc_tol))
  opts.srm.tol=1e-5; 
else
  opts.srm.tol=trunc_tol; 
end
if ( nargin<5 || isempty(max_order))
  opts.srm.max_ord=size(A,1); 
else
  opts.srm.max_ord=max_order; 
end

% some control settings for the LRADI 
if nargin<8 
    % ADI tolerance and maximum iteration number
    opts.adi.maxiter=100;
    opts.adi.res_tol=min( 1e-9, opts.srm.tol/100 );
    opts.adi.rel_diff_tol=1e-16;
    opts.norm='fro';
else
    if not(isfield(opts.adi,'maxiter')),  opts.adi.maxiter=100; end
    if not(isfield(opts.adi,'res_tol')), opts.adi.res_tol=min( 1e-9, opts.srm.tol/100 ); end
    if not(isfield(opts.adi,'rel_diff_tol')), opts.adi.rel_diff_tol=1e-16; end
    if not(isfield(opts,'norm')), opts.norm ='fro'; end
end 

% operations are done by the default set of user supplied finctions
oper = operatormanager('default');

%%
% Problem data

if not(issparse(E))||not(issparse(A))
    error('MESS:data', 'Both E and A need to be sparse.');
end
if sprank(E) < size(E,1)
    error('MESS:data', 'Only systems with invertible E are supported at the moment')
end

eqn.E_= E;
eqn.A_= A;
eqn.B = B;
eqn.C = C;
eqn.D = [];

n=oper.size(eqn, opts);

% Let us avoid E if it is actually the identity.
if norm(E-speye(n),'inf')==0
    eqn.haveE=0;
else
    eqn.haveE=1;
end


%%
% If not set outside, we use projection shifts
if not(isfield(opts.shifts,'method'))
    opts.shifts.num_desired=max(5, min(size(eqn.B,2), size(eqn.C,1)));
    opts.shifts.b0=ones(n,1);
    opts.shifts.method = 'projection';
end

opts.shifts.p=mess_para(eqn,opts,oper);

if opts.shifts.info && not(strcmp(opts.shifts.info, 'projection')), disp(opts.shifts.p); end
%%
% controllability
eqn.type='N';
outB = mess_lradi(eqn, opts, oper);
if outB.niter==opts.adi.maxiter
  warning('MESS:BT',['ADI did not converge for observability Gramian ' ...
                     'factor. Reduction results may be ' ...
                     'inaccurate']);
end

if info>1
  figure(1)
  semilogy(outB.res);
  title('AX + XA^T = -BB^T');
  xlabel('number of iterations');
  ylabel('normalized residual norm');
  drawnow
end

if info>0 
  disp('size outB.Z:')
  size(outB.Z)
end

%%
% observability

eqn.type = 'T';
outC = mess_lradi(eqn, opts, oper);
if outC.niter==opts.adi.maxiter
  warning('MESS:BT',['ADI did not converge for controllability Gramian ' ...
                     'factor. Reduction results may be ' ...
                     'inaccurate']);
end

if info>1
  figure(2)
  semilogy(outC.res);
  title('A^TX + XA = -C^TC');
  xlabel('number of iterations');
  ylabel('normalized residual norm');
  drawnow
end
if info>0
  disp('size outC.Z:')
  size(outC.Z)
end

%% execute square root method
[TL,TR,hsv] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);
%% compute ROM matrices
Ar = TL'*oper.mul_A(eqn, opts, 'N', TR, 'N');
Br = TL'*eqn.B;
Cr = eqn.C*TR;
Er = eye(size(Ar,1));
if isfield(eqn,'D')
    Dr=eqn.D;
else
    Dr = [];
end

%% if desired, plot the approximation results
if info>2
  if nargin<8
      opts.sigma.fmin = 1e-6;
      opts.sigma.fmax = 1e6;
      opts.sigma.nsample = 100;
  else
      if not(isfield(opts.sigma,'fmin')), opts.sigma.fmin = 1e-6; end
      if not(isfield(opts.sigma,'fmax')), opts.sigma.fmax = 1e6; end
      if not(isfield(opts.sigma,'nsample')), opts.sigma.nsample = 100; end
  end
  ROM = struct('A',Ar,'E',Er,'B',Br,'C',Cr,'D',Dr);
  mess_sigma_plot(eqn, opts, oper, ROM);
end

%% construct output information
if nargout > 4
    r  = size(Ar, 1);
    nr = size(A, 1) - length(hsv);
    
    outinfo = struct( ...
        'TL'      , TL, ...
        'TR'      , TR, ...
        'errbound', 2 * sum(hsv(r+1:end)) + nr * hsv(end),...
        'hsv'     , hsv);
end
