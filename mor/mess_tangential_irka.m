function [Er,Ar,Br,Cr,S,b,c,V,W] = mess_tangential_irka(varargin)
% The tangential IRKA method with automatic selection of initial shifts and
% tangential directions.
%
% Call
%   [Er,Ar,Br,Cr,S,b,c,V,W] = mess_tangential_irka(E,A,B,C,opts)
%
%   [Er,Ar,Br,Cr,S,c,V,W] = mess_tangential_irka(M,E,K,B,Cp,Cv,opts)
%
%   [Er,Ar,Br,Cr,S,c,V,W] = mess_tangential_irka(eqn,opts,oper)
%
% Inputs:
%  E,A,B,C    The mass, system, input and output matrices describing the
%             original system
%
%  M,E,K,B,   The mass, system, input and output matrices describing the
%  Cp,Cv      original system
%
%  opts       optional options structure with substructure 'irka'
%             
% Input fields in struct opts.irka:
%
%  r          The reduced order (optional, default: 25)
%
%  maxiter    maximum iteration number for the IRKA iteration 
%             (optional, default: 25)
%
%  shift_tol  bound for the relative change of the IRKA shifts used as
%             stopping criterion (optional, default: 1e-2)
%
%  h2_tol     bound for the relative change of the H2-norm compared to the
%             last stable ROM (optional, default: 100*eps)
%
%  info       0  : silent (default)
%             1  : print status info in each IRKA step
%             >1 : compute and show the sigma and error plots 
%                  (only for matrix inputs)
%
%  init       shift and direction initialization choice: (optional)
%              'subspace' chooses a random subspace and uses it to compute
%                         projected shifts and directions from the projected 
%                         EVP just like in the irka iteration.  (default)
%              'logspace' picks logspaced shifts in [0,1] and all ones
%                         as tangential directions
%              'random'   picks normally distributed random shifts and
%                         tangential directions
%              'rom'      an asymptotically stable initial guess for the
%                         reduced model of order r is given in opts.irka.Er,
%                         opts.irka.Ar, opts.irka.Br, opts.irka.Cr.
%
%  flipeig    flag marking whether or not to correct the signs of shifts in
%             the wrong halfplane. (optional, default: 1)
%
% Outputs:
%  Er,Ar,Br,Cr   The reduced system matrices.
%  S,b,c         The final shifts and tangential directions
%  V,W           The final projection matrices
%
% NOTE: Currently only standard state space systems and descriptor systems
% with E invertible are supported.
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


%% Choose usfs set if matrices were passed in
if nargin==5 % default first order system
    oper = operatormanager('default');
elseif nargin==7 % second order system
    oper = operatormanager('so_2');
end

%% Fill equation structure if matrices were passed in
if nargin==5
    eqn.A_ = varargin{2};
    if isempty(varargin{1})
        eqn.E_ = speye(size(varargin{2},1));
    else
        eqn.E_ = varargin{1};
        eqn.haveE = 1;
    end
    eqn.B = full(varargin{3});
    eqn.C = full(varargin{4});
    opts = varargin{5};
elseif nargin==7
    eqn.M_ = varargin{1};
    eqn.E_ = varargin{2};
    eqn.K_ = varargin{3};
    eqn.haveE = 1;
    eqn.B = [full(varargin{4}); zeros(size(varargin{4}))];
    eqn.C = [full(varargin{5}), full(varargin{6})];
    opts = varargin{7};
end

%%
if nargin==3
    eqn = varargin{1};
    opts = varargin{2};
    oper = varargin{3};
end

%% check field opts.irka
if not(isfield(opts,'irka')) || not(isstruct(opts.irka))
    opts.ikra = [];
end

if not(isfield(opts.irka,'r')) || isempty(opts.irka.r)
    opts.irka.r = 25;
end
if not(isfield(opts.irka,'maxiter')) || isempty(opts.irka.maxiter)
    opts.irka.maxiter = 25;
end
if not(isfield(opts.irka,'init')) || isempty(opts.irka.init)
    opts.irka.init = 'subspace';
end
if not(isfield(opts.irka,'shift_tol')) || isempty(opts.irka.shift_tol)
    opts.irka.shift_tol = 1e-2;
end
if not(isfield(opts.irka,'h2_tol')) || isempty(opts.irka.h2_tol)
    opts.irka.h2_tol = 100*eps;
end
if not(isfield(opts.irka,'info')) || isempty(opts.irka.info)
    opts.irka.info = 1;
end
if not(isfield(opts.irka,'flipeig')) || isempty(opts.irka.flipeig)
    opts.irka.flipeig = 1;
end

%% Initialize used usfs
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');
if not(result)
    error('MESS:data', 'Equation  data seems to be incomplete');
end
[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);

%% Initialization
n = oper.size(eqn, opts);
m = size(eqn.B,2);
p = size(eqn.C,1);

% new field saving whether the reduced order model is stable
isstab = zeros(opts.irka.maxiter,1);

% Choose the initial shifts and tangential dirctions
initial_rom=0;
switch opts.irka.init
    case 'subspace'
        U = get_initial_subspace(n,opts.irka.r);
        [T,S] = eig(full(U'*(oper.mul_A(eqn, opts, 'N', U, 'N'))),...
        full(U'*(oper.mul_E(eqn, opts, 'N', U, 'N'))));
        [ S, perm ] = mess_make_proper(diag(S));
        S = S';
        T = T(:,perm);
        b = (T\U'*eqn.B).';
        c = eqn.C*U*T;
    case 'random'
        S = abs(randn(opts.irka.r,1));
        b = randn(m,opts.irka.r);
        c = randn(p,opts.irka.r);
    case 'logspace'
        S = logspace(0,1,opts.irka.r);
        b = ones(m,opts.irka.r);
        c = ones(p,opts.irka.r);
    case 'rom'
        Er = opts.irka.Er;
        Ar = opts.irka.Ar;
        Br = opts.irka.Br;
        Cr = opts.irka.Cr;
        [T,S] = eig(Ar,Er);
        if any(real(diag(S))>=0)
            error(['The initial guess for the reduced system must ' ...
                   'be asymptotically stable!']); 
        end
        [ S, perm ] = mess_make_proper(diag(S));
        S = S';
        T = T(:,perm);
        b = (T\Br).';
        c = Cr*T;
        initial_rom=1;        
end
%% Start iteration
for iter = 1:opts.irka.maxiter
    %% save previous shifts
    S_old = S;
    %% save previous ROM
    if initial_rom || (iter>1 && isstab(iter-1))
      Erold = Er;
      Arold = Ar;
      Brold = Br;
      Crold = Cr;
    end
    %% Compute projection subspaces
    V = zeros(n,opts.irka.r);
    W = zeros(n,opts.irka.r);
    j = 1;
    while(j < opts.irka.r+1)        
        x = oper.sol_ApE(eqn, opts,'N',S(j),'N',eqn.B*b(:,j),'N');
        y = oper.sol_ApE(eqn, opts,'T',S(j),'T',eqn.C'*c(:,j),'N');
        if(imag(S(j)) ~= 0)
            V(:,j:j+1) = [real(x) imag(x)];
            W(:,j:j+1) = [real(y) imag(y)];
            j = j + 2;
        else
            V(:,j) = real(x);
            W(:,j) = real(y);
            j = j + 1;
        end
    end
    % find orthonormal bases
    [V,~] = qr(V,0);
    [W,~] = qr(W,0);
   
    %% Biorthogonalize V,W in the E inner product
    Er = (W'*oper.mul_E(eqn, opts, 'N', V, 'N'));
    [U,Sigma,Q] = svd(Er);
    Sigma=diag(ones(opts.irka.r,1)./sqrt(diag(Sigma)));
    W = W*U*Sigma;
    V = V*Q*Sigma;
    %% Compute new ROM

    Ar = W'*oper.mul_A(eqn, opts, 'N', V, 'N');
    Br = W'*eqn.B;
    Cr = eqn.C*V;    
    Er = eye(opts.irka.r); %by construction   
    %% Update interpolation points/tangential directions
    [T,S] = eig(Ar);
    [S, perm] = mess_make_proper(diag(S));
    S = S';
    T = T(:,perm);

    % make sure all shifts are in the coorect half plane.
    wrongsign=find(real(S)>0);
    if not(isempty(wrongsign))
        if (opts.irka.flipeig)
            warning('MESS:IRKA:unstable',...
            ['IRKA step %d : %d non-stable reduced eigenvalues have ' ...
             'been flipped.\n'], iter, length(wrongsign)); 
        else
            warning('MESS:IRKA:unstable',...
            'IRKA step %d : %d non-stable reduced eigenvalues detected.\n', ...
            iter, length(wrongsign));
        end
    else
        isstab(iter) = 1;
    end
    if (opts.irka.flipeig)
        S(wrongsign) = -S(wrongsign);
    end
    % update tangential directions
    b = (T\Br).';
    c = Cr*T;
    
    %% compute convergence indices
    % maximum pointwise relative change of the shifts
    shiftchg = norm(sort(S)-sort(S_old),'inf');
    % relative H2-change
    % Computation of the change only makes sense if we have an old model 
    % and both, the current and old, models are stable
    if (initial_rom && isstab(iter)) || ...
            ((iter>1) && (isstab(iter)&&any(isstab(1:iter-1)))) 
      romchg = mess_h2_rom_change(Er,Ar,Br,Cr,Erold,Arold,Brold,Crold,1);
    else
        if iter == 1
            romchg = 1.0;
        else
            romchg = Inf;
        end
    end
    %% If dsired print status message
    if opts.irka.info
        fprintf(['IRKA step %3d, rel. chg. shifts = %e , rel. H2-norm ' ...
                 'chg. ROM = %e\n'], iter, shiftchg, romchg); 
    end
    
    %% evaluate stopping criteria
    if(shiftchg < opts.irka.shift_tol) || (romchg < opts.irka.h2_tol)       
        break
    end
end
if ((iter == opts.irka.maxiter) && not((shiftchg < opts.irka.shift_tol)))
  warning('MESS:IRKA:convergence',...
      'IRKA: No convergence in %d iterations.\n', opts.irka.maxiter);
end

if opts.irka.info>1 && nargin>3
    ROM = struct('A',Ar,'E',Er,'B',Br,'C',Cr,'D',[]);
    if not(isfield(opts,'sigma')), opts.sigma = struct(); end
    if not(isfield(opts.sigma,'fmin')), opts.sigma.fmin=-6; end
    if not(isfield(opts.sigma,'fmax')), opts.sigma.fmax=6; end
    if not(isfield(opts.sigma,'nsample')), opts.sigma.nsample=100; end
    if not(isfield(opts.sigma,'info')), opts.sigma.info=opts.irka.info; end
    [~, eqn, opts, oper] = mess_sigma_plot(eqn,opts,oper,ROM);
end

[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
oper.sol_ApE_post(eqn, opts, oper);

