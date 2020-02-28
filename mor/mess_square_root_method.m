function [TL,TR,hsv,eqn,opts,oper] = mess_square_root_method(eqn,opts,oper,ZB,ZC)
% Square root method for the computation of the balanced and reduced
% system
%  
% Call
%  [TL,TR,hsv,eqn,opts,oper] = mess_square_root_method(eqn,opts,oper,ZB,ZC);
%
% Inputs:
%  eqn, opt, oper   the standard mess structures
%                   opts needs to have opts.srm.max_ord and opts.srm.tol
%                   set.
%  ZB, ZC           the (tall and skinny) Gramian factors
%
% Outputs: 
%  TL,TR            left and right truncation matrices
%  hsv              computed Hankel singular values
%
% The implementation (especially in the case E!=I) follows the
% derivation in:
%	Efficient Numerical Solution of Large Scale Algebraic Matrix
%	Equations in PDE Control and Model Order Reduction; 
%   Saak, Jens;
%   Dissertation, TU Chemnitz; 2009. 
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
%% Check necessary control data
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');
[eqn, opts, oper] = oper.mul_E_pre(eqn,opts,oper);
if not(result)
    error('MESS:control_data',...
        'system data is not completely defined or corrupted');
end
if isfield(opts,'srm')
   if not(isfield(opts.srm,'tol'))
       error('Missing truncation tolerance opts.srm.tol');
   end
   if not(isfield(opts.srm,'max_ord'))
        opts.srm.max_ord=oper.size(eqn,opts);
   end
   if not(isfield(opts.srm,'info'))
       opts.srm.info = 0;
   end
else
    error('Missing srm substructure in opts argument.');
end

%% Compute SVD of Gramian factor product in the correct inner product
if eqn.haveE
    [U0,S0,V0] = svd(ZC'*oper.mul_E(eqn, opts, 'N', ZB, 'N'),0);
else
    [U0,S0,V0] = svd(ZC'*ZB,0);
end

%% Determine possible and desired ROM order
hsv=diag(S0);
ks=length(hsv);
nr=oper.size(eqn,opts)-ks;
k=ks;
while (2*sum(hsv(k-1:ks))+nr*hsv(ks)<opts.srm.tol)&&(k>2)
    k=k-1;
end
k0=k;

r= min([opts.srm.max_ord k0]);
if opts.srm.info>0
    fprintf(1,['reduced system order: %d',...
        '  (max possible/allowed: %d/%d)\n\n'],r,ks,opts.srm.max_ord);
end

% Compute the truncating projection matrices

S = sparse(1:r, 1:r, 1 ./ sqrt(hsv(1:r)));

TL = (ZC*U0(:,1:r))*S;
TR = (ZB*V0(:,1:r))*S;

% augment projection matrices by preselected columns
if isfield(opts.srm,'V') && ismatrix(opts.srm.V)
    TL = [TL,opts.srm.V];
end
if isfield(opts.srm,'W') && ismatrix(opts.srm.W)
    TR = [TR, opts.srm.W];
end

[eqn, opts, oper] = oper.mul_E_post(eqn,opts,oper);
