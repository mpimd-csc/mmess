function y=riccati(Z,x,eqn,oper,opts,D)
% Computes matrix vector product with the Riccati operator.
%  y = A^T*ZZ^T*x + ZZ^T*Ax + C^T*Cx - ZZ^T*BB^T*ZZ^T*x  or
%  y = A^T*ZZ^T*Ex + C^T*Cx + E^T*ZZ^T*Ax - BB^T*ZZ^Tx
%
% Input:
%  Z         Low-rank solution factor of the Riccati equation
%  x         vector for matrix vector product
%  eqn       structure with data for A, E and fields B, C
%                  eqn.E(optional, eqn.haveE specifies whether it is
%                  there) in the above equation with ZZ' approximating X
%  oper      structure contains function handles for operations with
%                  A, E
%  opts    structure contains parameters for the algorithm
%
%  D       solution factor D for LDL^T formulation in case opts.LDL_T=1
%
% Output:
%  y        result of matrix vector product
%
% generalized equations
% y = A'*Z*Z'*E*x + C'*C*x + E'*Z*Z'*(A*x - B*B'*Z*Z'*x)
% y = A*Z*Z'*E'*x + B*B'*x + E*Z*Z'*(A'*x - C'*C*Z*Z'*x)
%
% or
%
% y = A'*Z*D*Z'*E*x + C'*eqn.S*C*x
%   + E'*(Z*D*Z'*A*x - B*B'*Z*D*Z'*x)
%
% uses operatorfunctions mul_E, mul_A, mul_ApE (inside BDF methods)

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


%% Check input
if nargin < 6 && opts.LDL_T
    error('for the LDL^T version the information to build D must be passed');
end

if not(isfield(opts,'bdf')), opts.bdf=[]; end
if isstruct(opts.bdf) && isfield(opts.bdf, 'tau') && isfield(opts.bdf, 'beta')
    bdf = 1;
    pc = -1 / (2 * opts.bdf.tau * opts.bdf.beta);
else
    bdf = 0;
end

if isempty(D) && opts.LDL_T
    error('LDL^T formulation needs D to get passed.');
end

%% Compute MVP
% compute z = X*E*x, z=X*E'*x, or simply z = X*x
if eqn.haveE
    if opts.LDL_T
        if eqn.type == 'T'
            z=Z*mess_LDL_mul_D(eqn, D, Z'*(oper.mul_E(eqn, opts,'N',x,'N')));
        else
            z=Z*mess_LDL_mul_D(eqn, D, Z'*(oper.mul_E(eqn, opts,'T',x,'N')));
        end
    else
        if eqn.type == 'T'
            z=Z*(Z'*(oper.mul_E(eqn, opts,'N',x,'N')));
        else
            z=Z*(Z'*(oper.mul_E(eqn, opts,'T',x,'N')));
        end
    end
else
    if opts.LDL_T
        z=Z*mess_LDL_mul_D(eqn, D, Z'*x);
    else
        z=Z*(Z'*x);
    end
end

% Compute y1 = A'*z + C'*(C*x)  or y1 =  A*z + B*B'*x  
% and     y2 = X*(A*x - B*B'*z) or y2 =  X*(A'*x - C'*C*z)
if bdf
    if eqn.type == 'T'
        y1 = (opts.bdf.tau * opts.bdf.beta) * ...
            oper.mul_ApE(eqn, opts,'T',pc, 'T', z,'N') + eqn.pC'*eqn.S*(eqn.pC*x);
        y2 = Z*mess_LDL_mul_D(eqn, D, Z'*((opts.bdf.tau * opts.bdf.beta) * ...
            oper.mul_ApE(eqn, opts,'N',pc, 'N', x,'N') ...
            - (opts.bdf.tau * opts.bdf.beta) * eqn.B*(eqn.B'*z)));
    else
        y1 = (opts.bdf.tau * opts.bdf.beta) * ...
            oper.mul_ApE(eqn, opts,'N',pc, 'N', z,'N') + eqn.pB*eqn.S*(eqn.pB'*x);
        y2 = Z*mess_LDL_mul_D(eqn, D, Z'*((opts.bdf.tau * opts.bdf.beta) * ...
            oper.mul_ApE(eqn, opts,'T',pc, 'T', x,'N')...
            - (opts.bdf.tau * opts.bdf.beta) * eqn.C'*(eqn.C*z)));
    end
elseif eqn.haveUV && eqn.sizeUV1
    if eqn.type == 'T'
        y1 = oper.mul_A( eqn, opts, 'T', z, 'N' ) ...
            + eqn.pV(:,1:eqn.sizeUV1) * (eqn.U(:,1:eqn.sizeUV1)' * z) ...
            + eqn.pC' * (eqn.pC * x);
        y2 = Z * (Z'*(oper.mul_A( eqn, opts, 'N', x, 'N' ) ...
            + eqn.U(:,1:eqn.sizeUV1) * (eqn.V(:,1:eqn.sizeUV1)' * x) ...
            - eqn.B * (eqn.B' * z)));
    else
        y1 = oper.mul_A( eqn, opts, 'N', z, 'N' ) ...
            + eqn.pU(:,1:eqn.sizeUV1) * (eqn.V(:,1:eqn.sizeUV1)' * z) ...
            + eqn.pB * (eqn.pB' * x);
        y2 = Z * (Z'*(oper.mul_A( eqn, opts, 'T', x, 'N' ) ...
            + eqn.V(:,1:eqn.sizeUV1) * (eqn.U(:,1:eqn.sizeUV1)' * x) ...
            - eqn.C' * (eqn.C * z)));
    end
else
    if opts.LDL_T
        if eqn.type == 'T'
            y1 = oper.mul_A(eqn, opts,'T',z,'N') + eqn.pC'*eqn.S*(eqn.pC*x);
            y2 = Z*mess_LDL_mul_D(eqn, D, Z'*(oper.mul_A(eqn, opts,'N',x,'N')-eqn.B*(eqn.B'*z)));
        else
            y1 = oper.mul_A(eqn, opts,'N',z,'N') + eqn.pB*eqn.S*(eqn.pB'*x);
            y2 = Z*mess_LDL_mul_D(eqn, D, Z'*(oper.mul_A(eqn, opts,'T',x,'N')-eqn.C'*(eqn.C*z)));
        end
    else
        if eqn.type == 'T'
            y1 = oper.mul_A(eqn, opts,'T',z,'N') + eqn.pC'*(eqn.pC*x);
            y2 = Z*(Z'*(oper.mul_A(eqn, opts,'N',x,'N')-eqn.B*(eqn.B'*z)));
        else
            y1 = oper.mul_A(eqn, opts,'N',z,'N') + eqn.pB*(eqn.pB'*x);
            y2 = Z*(Z'*(oper.mul_A(eqn, opts,'T',x,'N')-eqn.C'*(eqn.C*z)));
        end
    end
end

% Now y = (C'*C + A'*X*E+E'*X*A - E'*X*B*B'*X*E)*x,  
%     y = (B*B' + A*X*E'+E*X*A' - E'*X*C'*C*X*E)*x  
%     or the same with E=I
if eqn.haveE
    if eqn.type == 'T'
        y  = y1 + oper.mul_E(eqn, opts,'T',y2,'N');
    else
        y  = y1 + oper.mul_E(eqn, opts,'N',y2,'N');
    end
else
    y = y1 + y2;
end
