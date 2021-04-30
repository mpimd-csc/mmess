function [Q, R] = mess_mgs(A,E)
%
%  function [Q, R] = mess_mgs(A,E);
%
% modified Gram-Schmidt orthogonalization of the columns of
% A. The columns of A are assumed to be linearly independent.
%
% [Q, R] = mgs(A) returns a matrix Q with orthonormal columns
% and an invertible upper triangular matrix R so that A = Q*R.
% If only orthogonalization is wanted, the accumulation of R can be
% avoided by ommiting the output parameter:
% Q = mgs(A)
%
% [Q, R] = mgs(A,E) performes the orthonomalization in the E scalar
% product, i.e., E needs to be symmetric positive definite.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isnumeric(A)) || not(ismatrix(A))
    error('MESS:error_arguments','A has to be a matrix')
end
if (nargin == 2) && (not(isnumeric(E)) || not(ismatrix(E)))
    error('MESS:error_arguments','E has to be a matrix')
end
if (nargin == 2) && (size(A, 1) ~= size(E, 2))
    error('MESS:error_arguments','number of columns of E differs with number of rows of A');
end

  [m,n]=size(A);

  if (nargin==2)&&(any(any(E'-E)))
    error('MGS:input matrix E needs to be selfadjoint');
  end

  R = zeros(n,n);
  Q = zeros(m,n);
  for k=1:n
    if nargout==2
      if nargin==2
        R(k,k) = Enorm(E,A(:,k));
      else
        R(k,k) = norm(A(:,k));
      end
      Q(:,k) = A(:,k)/R(k,k);
      for j=k+1:n
        if nargin==2
          R(k,j) = Q(:,k)'*(E*A(:,j));
        else
          R(k,j) = Q(:,k)'*A(:,j);
        end
        A(:,j) = A(:,j)-Q(:,k)*R(k,j);
      end
    else
      if nargin==2
        R = Enorm(E,A(:,k));
      else
        R = norm(A(:,k));
      end
      Q(:,k) = A(:,k)/R;
      for j=k+1:n
        if nargin==2
          R = Q(:,k)'*(E*A(:,j));
        else
          R = Q(:,k)'*A(:,j);
        end
        A(:,j) = A(:,j)-Q(:,k)*R;
      end
    end
  end

  function nrm=Enorm(E,x)
    nrm=sqrt(x'*(E*x));
