function Y = mess_LDL_mul_D(eqn,D,X)
%
%  Computes Y = D * X
%
% Here D is the center matrix in the LDL' representation. It can be given 
% either as the full D matrix such that D * X can be computed directly, or
% implicitly stored, e.g., as a vector D such that the actual D matrix is
% kron(diag(D), eqn.S). The later is done, e.g, in the LDL' formulation of 
% the low-rank ADI, while in differnetial Riccati solvers we typically have
% a small D matrix such that the actual D is kron(D,eqn.S)
%
% In the latter cases we exploit 
%
%   kron(F,eqn.S)*vec(Z) = vec(eqn.S * Z * F')  
%
% where F is either diag(D) or D, which can be efficiently expressed using
% reshape. 

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

% Do a consistency check on the inputs (no full check is performed though)
ns = size( eqn.S_diag, 1 );

if isvector(D), D=diag(D); end
nd = size(D,1); 
[mx, nx] = size(X);

if mx~=nd && mx~=ns*nd
    error(['The number of rows in X must be either equal to the ' ...
           'number of columns in D or the product of the numbers ' ...
           'of columns in D and eqn.S']);  
end

if mx==ns*nd % the full D is actually kron(diag(D),eqn.S)    
    % turn each column in X into an ns by nd matrix
    X = reshape(X,ns,nd,nx);
    % allocate Y as a same size 3d array
    Y = zeros(size(X));
    for i=1:nx
        Y(:,:,i) =diag(eqn.S_diag)*(X(:,:,i)*D');
    end
    % turn the matrified result into columns again.
    Y = reshape(Y,mx,nx); 
else % D is the full D matrix in the LDL^T decomposition
    Y = D * X;
end