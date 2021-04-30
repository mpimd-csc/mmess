function [Er,Ar,Br,Cr,outinfo]=mess_balanced_truncation(varargin)
% Lyapunov Balanced truncation for descriptor systems with invertible E.
%
%  [Er,Ar,Br,Cr,outinfo] = mess_balanced_truncation(E,A,B,C,max_order,trunc_tol,info,opts)
%                               or
%  [Er,Ar,Br,Cr,outinfo] = mess_balanced_truncation(sys,max_order,trunc_tol,info,opts)
%                        with sys = sparss(A,B,C,D,E)
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
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%



%%
if isa(varargin{1}, 'sparss')
    if nargin> 1
    max_order = varargin{2};
    end

    if nargin > 2
    trunc_tol = varargin{3};
    end

    if nargin > 3
    info = varargin{4};
    end

    if nargin > 4
    opts = varargin{5};
    end
else
    if nargin > 4
        max_order = varargin{5};
    end

    if nargin > 5
        trunc_tol = varargin{6};
    end

    if nargin > 6
         info = varargin{7};
    end

    if nargin > 7
         opts = varargin{8};
    end
end

% control verbosity of the computations
if nargin<7
  opts.adi.info=0;
  opts.bt.info=0;
  opts.srm.info=0;
  info=0;
  opts.shifts.info = 0;
else
    if nargin < 8
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


%% check if sparss or matrices were passed in
if isa(varargin{1}, 'sparss')
    [eqn, oper] = mess_wrap_sparss(varargin{1});
    n=oper.size(eqn, opts);
    if exist('eqn.D', 'var')
       warning('MESS:ignored',...
                'Matrix D is supposed to be empty. Data is ignored.');
    end
    % save non truncuated matrices
    B = eqn.B;
    C = eqn.C;
    % Note that we truncate B and C for the best robustness of the low-rank
    % Lyapunov solvers, here.
    eqn.B = mess_column_compression(eqn.B, 'N');
    eqn.C = mess_column_compression(eqn.C, 'T');
    eqn.D = [];
else
    % Problem data
    if not(issparse(varargin{1}))||not(issparse(varargin{2}))
        error('MESS:data', 'Both E and A need to be sparse.');
    end
    % needs to be saved in extra variables for octave 4.2 compability
    rank_1 = sprank(varargin{1});
    dimension_1 = size(varargin{1},1);
    if rank_1 < dimension_1
        error('MESS:data', 'Only systems with invertible E are supported at the moment')
    end

    eqn.E_= varargin{1};
    eqn.A_= varargin{2};
    % save non truncuated matrices
    B = varargin{3};
    C = varargin{4};
    % Note that we truncate B and C for the best robustness of the low-rank
    % Lyapunov solvers, here.
    eqn.B = mess_column_compression(B, 'N');
    eqn.C = mess_column_compression(C, 'T');
    eqn.D = [];

    % operations are done by the default set of user supplied finctions
    oper = operatormanager('default');
    n=oper.size(eqn, opts);

    % Let us avoid E if it is actually the identity.
    if norm(eqn.E_-speye(n),'inf')==0
        eqn.haveE=0;
    else
        eqn.haveE=1;
    end
end

% BT tolerance and maximum order for the ROM
if (nargin < 6 || isempty(trunc_tol))
  opts.srm.tol=1e-5;
else
  opts.srm.tol=trunc_tol;
end
if ( nargin < 5 || isempty(max_order))
  opts.srm.max_ord=size(eqn.A_,1);
else
  opts.srm.max_ord=max_order;
end

% some control settings for the LRADI
if nargin < 8
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


%%
% If not set outside, we use projection shifts
if not(isfield(opts.shifts,'method'))
    opts.shifts.num_desired=max(5, min(size(eqn.B,2), size(eqn.C,1)));
    opts.shifts.b0=ones(n,1);
    opts.shifts.method = 'projection';
end

opts.shifts.p = mess_para(eqn, opts, oper);

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

if info > 1
  figure(1)
  semilogy(outB.res);
  title('AX + XA^T = -BB^T');
  xlabel('number of iterations');
  ylabel('normalized residual norm');
  drawnow
end

if info > 0
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

if info > 1
  figure(2)
  semilogy(outC.res);
  title('A^TX + XA = -C^TC');
  xlabel('number of iterations');
  ylabel('normalized residual norm');
  drawnow
end
if info > 0
  disp('size outC.Z:')
  size(outC.Z)
end

%% execute square root method
[TL,TR,hsv] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);
%% compute ROM matrices
% Note that we use the original B and C since the ones in eqn have been
% truncated above.
Ar = TL'*oper.mul_A(eqn, opts, 'N', TR, 'N');
Br = TL'*B;
Cr = C*TR;
Er = eye(size(Ar,1));

if isfield(eqn,'D')
    Dr=eqn.D;
else
    Dr = [];
end

%% if desired, plot the approximation results
if info > 2
  if nargin < 8
      opts.sigma.fmin = 1e-6;
      opts.sigma.fmax = 1e6;
      opts.sigma.nsample = 100;
  else
      if not(isfield(opts.sigma,'fmin')), opts.sigma.fmin = 1e-6; end
      if not(isfield(opts.sigma,'fmax')), opts.sigma.fmax = 1e6; end
      if not(isfield(opts.sigma,'nsample')), opts.sigma.nsample = 100; end
  end
  ROM = struct('A',Ar,'E',Er,'B',Br,'C',Cr,'D',Dr);
  % for the evaluations here we need the original B and C in the eqn
  % structure.
  eqn.B = B;
  eqn.C = C;
  mess_sigma_plot(eqn, opts, oper, ROM);
end

%% construct output information
if nargout > 4
    r  = size(Ar, 1);
    nr = size(eqn.A_, 1) - length(hsv);

    outinfo = struct( ...
        'TL'      , TL, ...
        'TR'      , TR, ...
        'errbound', 2 * sum(hsv(r+1:end)) + nr * hsv(end),...
        'hsv'     , hsv);
end
