function [nrm,k,T,V,eqn,fopts,oper] = mess_res2_norms(Z,Rmul,eqn,fopts,oper,resopts,D)
% Computes the 2 Norm of the residual of Z for the symmetric operator given 
%   by y=Rmul(Z, x, eqn, opts,oper);
% e.g., for the generalized Lyapunov residual Rmul implements
%   R := FZZ^T*E^T + EZZ^T*F^T + GG^T (1) 
% 
%
% That means res2 computes the spectral radius of this operator R by a
% Lanzcos iteration. The function Rmul should exploite the structure of F 
% and rectangular structure of Z and G. Thus it can be computed in O(n) 
% effort and is therefore much cheaper than the computation of the e.g. the
% Frobenius norm.
%
% Method:
%
%    The Lanczos method produces matrices V and T such that
%
%      V(:,1) in span{rv},
%      V'*V = eye(k+1),
%      F*V(:,1:k) = V*T.
%
%  Remark:
%    V is only constructed if its required by the output arguments, or if
%    reorthogonalization is used.
% 
%    This implementation does not check for (near-)breakdown!
%
%  
% Input:                 
%  Z                Low-rank solution factor of operator
%
%  eqn              structure with data for operator 
%
%  Rmul             function handle to a function implementing the 
%                   multiplication with the residual operator of interest.                
%
%  fopts            full options structure (passed on to function handles in
%                   oper)
%
%  oper             structure contains function handles for operations with
%                   data in eqn
%
%  resopts          options structure with fields
%  resopts.maxiter  maximal number of Arnoldi steps (usually k<<n)
%                   (optional - chosen as 10 if omitted)
%  resopts.tol      relative accuracy of the norm computation       
%                   (optional - chosen as 1e-6 if omitted)
%  resopts.rv       initial n-vector
%                   (optional - chosen by random, if omitted)
%  resopts.orth     reorthogonalization flag
%                   (optional - switched off, if omitted)
%  D                solution factor D for LDL^T formulation in case
%                   fopts.LDL_T = 1 
%
%
% Output:
%  nrm              the residual of the iterate X=Z*Z' or X = Z*D*Z'
%
%  k                number of Lanzcos steps taken
%
%  T                matrix T ((k+1)-x-k matrix, symmetric tridiagonal);
%
%  V                matrix V (n-x-(k+1) matrix, orthogonal columns).
%
% uses eventually operatorfunctions in Rmul

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

%% check Rmul
[~, eqn, fopts, oper] = oper.init(eqn,fopts,oper,'A','E');
[eqn, fopts, oper] = oper.init_res_pre(eqn,fopts,oper);
[eqn, fopts, oper] = oper.mul_A_pre(eqn,fopts,oper);
[eqn, fopts, oper] = oper.mul_E_pre(eqn,fopts,oper);

n = size(Z,1); % Get system order.
switch Rmul
    case 'lyapunov'
        Rmul=@lyapunov;
    case 'riccati'
        Rmul=@riccati;
    otherwise
        error('MESS:control_data','Rmul has to be ''lyapunov'' or ''riccati''');
end

if not(isfield(resopts,'res'))
  error('MESS:control_data','residual control structure res missing.');
else
    res=resopts.res;
end

if not(isfield(res,'maxiter'))||isempty(res.maxiter)
    warning('MESS:control_data','res.maxiter is set to 10 (default)');
    res.maxiter = 10;
else
    if res.maxiter >= n-1
        error('maxiter must be smaller than the order of A!');
    end
end

if not(isfield(res,'tol'))||isempty(res.tol)
    warning('MESS:control_data','res.tol is set to 1e-6 (default)');
    res.tol= 1e-6; 
end

if not(isfield(res,'rv'))||isempty(res.rv)
    res.rv = oper.init_res(eqn, fopts, oper, randn(n,1)); 
else
    res.rv =  oper.init_res(eqn, fopts, oper, res.rv);
end

if not(isfield(res,'orth'))||isempty(res.orth)
    warning('MESS:control_data','res.orth is set to 0 (default)');
    res.orth=0; 
end

if not(isfield(fopts,'LDL_T')), fopts.LDL_T = 0; end
if not(isfield(eqn,'haveUV')), eqn.haveUV = 0; end

if nargin< 7, D=[]; end

created_projected_data = 0;
if eqn.type=='N'
    if not(isfield(eqn,'pB')) || (eqn.haveUV && not(isfield(eqn, 'pU')))
        created_projected_data = 1;
        if eqn.haveUV
            [W, ~, eqn, fopts, oper] = oper.init_res(eqn, fopts, oper, [eqn.B, eqn.U]);
            eqn.pB = W(:,1:size(eqn.B,2));
            eqn.pU = W(:,size(eqn.B,2)+1:end);
        else
            [pBtemp, ~, eqn, fopts, oper] = oper.init_res(eqn, fopts, oper, eqn.B);
            eqn.pB = pBtemp;
        end
    end
else
    if not(isfield(eqn,'pC')) || (eqn.haveUV && not(isfield(eqn, 'pV')))
        created_projected_data = 1;
        if eqn.haveUV
            [W, ~, eqn, fopts, oper] = oper.init_res(eqn, fopts, oper, [eqn.C', eqn.V]);
            eqn.pC = W(:,1:size(eqn.C,1))';
            eqn.pV = W(:,size(eqn.C,1)+1:end);
        else
            [pCtemp, ~, eqn, fopts, oper] = oper.init_res(eqn, fopts, oper, eqn.C');
            eqn.pC = pCtemp';
        end
    end
end

%% start computation
T = zeros(res.maxiter+1,res.maxiter);
v1 =(1.0/norm(res.rv))*res.rv;

% initial step
% Matrix-vector product R*v1
w=Rmul(Z,v1,eqn,oper,fopts,D);

T(1,1)=v1'*w;
r=w-T(1,1)*v1;
T(1,2)=norm(r);
v2 = r./T(1,2);
T(2,1)=T(1,2)';

% store vectors in V if needed
if res.orth || nargout>3
    V = zeros(n,res.maxiter+1);
    V(:,1) = v1;
    V(:,2)=v2;
end

nrm = 0;

for k=2:res.maxiter-1
    %Lanczos 3-term recursion 
    
    % Matrix-vector product R*v2
    w = Rmul(Z,v2,eqn,oper,fopts,D);

    T(k,k)=v2'*w;
    r=w-T(k-1,k)*v1-T(k,k)*v2;
    
    %re-orthogonalization by MGS
    if res.orth
        for j=1:k
            r = r - (V(:,j)'*r)*V(:,j);
        end
    end
    
    T(k,k+1)=norm(r);
    v1=v2;
    v2 = r./T(k,k+1);
    if res.orth || nargout>3
        V(:,k+1)=v2;
    end
    T(k+1,k)=T(k,k+1)';
    nrmold= nrm;
    nrm = max(abs(eig(T(1:k,1:k))));
    if abs(nrm-nrmold)<res.tol*nrm
        break;
    end
end

if created_projected_data
    if eqn.type == 'N'
        if eqn.haveUV, eqn = rmfield(eqn,'pU'); end
        eqn = rmfield(eqn,'pB'); 
    else
        if eqn.haveUV, eqn = rmfield(eqn,'pV'); end
        eqn = rmfield(eqn,'pC'); 
    end
end
[eqn, fopts, oper] = oper.init_res_post(eqn,fopts,oper);
[eqn, fopts, oper] = oper.mul_A_post(eqn,fopts,oper);
[eqn, fopts, oper] = oper.mul_E_post(eqn,fopts,oper);
