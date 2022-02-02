function Y = mess_LDL_mul_D(eqn,D,X)
%
%  Computes Y = D * X
%
% Here D is the center matrix in the LDL' representation. It can be given
% either as the full D matrix such that D * X can be computed directly, or
% implicitly stored, e.g., as a vector D such that the actual D matrix is
% kron(diag(D), eqn.S). The later is done, e.g, in the LDL' formulation of
% the low-rank ADI, while in differential Riccati solvers we typically have
% a small D matrix such that the actual D is kron(D,eqn.S)
%
% In the latter cases we exploit
%
%   kron(F,eqn.S)*vec(Z) = vec(eqn.S * Z * F')
%
% where F is either diag(D) or D, which can be efficiently expressed using
% reshape.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


% Do a consistency check on the inputs (no full check is performed though)
ns = size( eqn.S_diag, 1 );

if isvector(D), D = diag(D); end

nd = size(D,1);

[mx, nx] = size(X);

if (mx ~= nd) && (mx ~= ns*nd) % TODO This could be an assert?
    error(['MESS:mess_LDL_mul_D:The number of rows in X must be either equal to ' ...
           'the number of columns in D or the product of the numbers ' ...
           'of columns in D and eqn.S']);
end

if mx == (ns*nd) % the full D is actually kron(diag(D),eqn.S)
    % turn each column in X into an ns by nd matrix
    X = reshape(X,ns,nd,nx);
    % allocate Y as a same size 3d array
    Y = zeros(size(X));
    for k=1:nx
        Y(:,:,k) =diag(eqn.S_diag)*(X(:,:,k)*D');
    end
    % turn the matrified result into columns again.
    Y = reshape(Y,mx,nx);
else % D is the full D matrix in the LDL^T decomposition
    Y = D * X;
end
