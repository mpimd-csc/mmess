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
% avoided by omitting the output parameter:
% Q = mgs(A)
%
% [Q, R] = mgs(A,E) performs the orthonormalization in the E scalar
% product, i.e., E needs to be symmetric positive definite.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isnumeric(A)) || not(ismatrix(A))
    error('MESS:error_arguments','A has to be a matrix');
end

if (nargin == 2)

    if not(isnumeric(E)) || not(ismatrix(E))
        error('MESS:error_arguments','E has to be a matrix');
    end

    if size(A, 1) ~= size(E, 2)
        error('MESS:error_arguments','number of columns of E differs with number of rows of A');
    end

    if any(any(E'-E))
        error('MGS:input matrix E needs to be selfadjoint');
    end
else
    E = 1;
end

[m,n] = size(A);

R = zeros(n,n);
Q = zeros(m,n);

for k=1:n
    if nargout == 2

        R(k,k) = Enorm(E,A(:,k));

        Q(:,k) = A(:,k)/R(k,k);

        for l=k+1:n

            R(k,l) = Edot(Q(:,k)',A(:,l));

            A(:,l) = A(:,l)-Q(:,k)*R(k,l);
        end
    else

        R = Enorm(E,A(:,k));

        Q(:,k) = A(:,k)/R;

        for l=k+1:n

            R = Edot(Q(:,k)',A(:,l));

            A(:,l) = A(:,l)-Q(:,k)*R;
        end
    end
end
end

function mm = Edot(E,x,y)

    if isscalar(E) && (E == 1)

        mm = x * y;
    else

        mm = x * (E*y);
    end
end

function nrm = Enorm(E,x)

    if isscalar(E) && (E == 1)

        nrm = norm(x);
    else

        nrm = sqrt(x'*(E*x));
    end
end

