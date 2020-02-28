function y=lyapunov(Z,x,eqn,oper,opts, D)
% Computes matrix vector product with the Lyapunov operator.
%
% Input:           
%  Z         Low-rank solution factor of the Riccati equation
%  x         vector for matrix vector product
%  eqn       structure with data for A, E and fields G  
%                  eqn.E(optional, eqn.haveE specifies whether it is
%                  there) in the above equation with ZZ' or LDL' approximating X
%                  eqn.haveUV specfies whether feedback is there
%                  
%  oper      structure contains function handles for operations with
%                  A, E
%  opts      full options structure (passed on to function handles in oper) 
%
%  D         solution factor D for the LDL^T
%            decomposition, i.e., opts.LDL_T=1
%  eqn.type 'N' or 'T' for type of Lyapunov equation
%
% Output:
%  y        result of matrix vector product

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

if eqn.type=='N'
  adjoint='T';
else
  adjoint='N';
end

if isempty(D) && opts.LDL_T
    error('LDL^T formulation needs D to get passed.'); 
end

if eqn.haveE
    if opts.LDL_T
        z=Z*mess_LDL_mul_D(eqn,D, Z'*oper.mul_E(eqn, opts,adjoint,x,'N'));
    else
        z=Z*(Z'*(oper.mul_E(eqn, opts,adjoint,x,'N')));
    end
else
    if opts.LDL_T
        z=Z*mess_LDL_mul_D(eqn, D, Z'*x);
    else
        z=Z*(Z'*x);
    end
end

y1 = oper.mul_A(eqn, opts, eqn.type, z, 'N');
y2 = oper.mul_A(eqn, opts, adjoint, x, 'N');

if eqn.haveUV
    if eqn.type=='N'
        y1 = y1 + eqn.pU*(eqn.V'*z);
        y2 = y2 + eqn.V*(eqn.U'*x);
    else
        y1 = y1 + eqn.pV*(eqn.U'*z);
        y2 = y2 + eqn.U*(eqn.V'*x);
    end
end
if opts.LDL_T
    y2=Z*mess_LDL_mul_D(eqn, D, Z'*y2);
else
    y2=Z*(Z'*y2);
end

if eqn.haveE
    y2= oper.mul_E(eqn, opts,eqn.type, y2, 'N');
end

if opts.LDL_T
    y=y1+y2+eqn.G*(eqn.S*(eqn.G'*x));
else
    y=y1+y2+eqn.G*(eqn.G'*x);
end

% in case of Rosenbrock we get a -1/(2*timestep)*(E'*Z*Z'*E)
% from both F and F'
if isfield(opts,'rosenbrock')&&not(isempty(opts.rosenbrock))
  if eqn.haveE       % generalized equations
    y=y-(1/opts.rosenbrock.stepsize)*oper.mul_E(eqn, opts,eqn.type,z,'N');
  else
    y=y-(1/opts.rosenbrock.stepsize)*z;
  end
end
