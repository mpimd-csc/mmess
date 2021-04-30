function [Y] = care_nwt_fac(Y0,A,B,C,tol,maxsteps)
% Newton's method for continuous-time algebraic Riccati equations
%  (CARE)  0  =  C'C  +  A' X  +  X A  -  X BB' X =: R(X)
%                                               X = Y' Y
%
%
% Input:
%  Y0        initial starting guess s.t. A - BB'Y_0'Y_0 is stable.
%                  Note: this is not checked - if Y_0 is not stabilizing, then
%                  the iteration may fail to converge or converge to a
%                  non-stabilizing solution!
%  A         Matrix from (CARE)
%  B         Matrix from (CARE)
%  C         Matrix from (CARE)
%  tol       stopping criterion, i.e. the iteration is stopped if
%                  |R(X)|_F/max(1,|Y'Y|_F) <= tol
%                  , default sqrt(eps*n)
%
%  maxsteps  maximum number of iteration steps, default 50
%
% Output:
%  Y        approximate factor solution of CARE so that X=Y'Y

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

narginchk(4,6);
n = size(A,1);

%% check matrix sizes
if size(A,2) ~= n
    error('A must be square.');
end
if size(B,1) ~= n
    error('B must have the same number of rows as A.');
end
if size(C,2) ~= n
    error('C must have the same number of columns as A.');
end

if nargin < 6,  maxsteps = 50;  end
if (nargin < 5)
  tol = sqrt(eps*n);
else
    if tol < sqrt(eps*n)
        tol = sqrt(eps*n);
        warning('MESS:CARE_NWT_FAC',...
                'Error tolerance too small, may not be achieved!');
    end
end

%% Initialization
iter = 0;
if isempty(Y0)
    Y = zeros(1,n);
else
    if size(Y0,2) ~= n
        error('Y0 must have the same number of rows as A.');
    else
        Y = Y0;
    end
end
YA    = Y*A;
YB    = Y*B;
nres  = norm(C'*C + YA'*Y + Y'*YA - YB*YB','fro');
Xnorm = norm(Y*Y','fro');
Err = nres/max(1,Xnorm);
onemore     = 0;
convergence = Err <= tol;

%% Newton iteration
while (iter < maxsteps) && ((not(convergence)) || (convergence && (onemore < 2)))
  % Here one may employ RRQR to compress W.
  W        = [C; YB'*Y];
  Y        = lyap_sgn_fac(A - B*(YB)'*Y,W);
  YA       = Y*A;
  YB       = Y*B;
  nres     = norm(C'*C + YA'*Y + Y'*YA - (Y'*YB)*(YB'*Y),'fro');
  Xnorm    = norm(Y*Y','fro');
  iter = iter + 1;
% Uncomment next line for verbose mode.
%  fprintf('||R(X_%i)||/||X|| = %d\n', iter, nres/Xnorm)
  Err = nres/max(1,Xnorm);
  convergence = Err <= tol;
  if convergence,  onemore = onemore + 1;  end
end

if (iter == maxsteps) && (nres/max(1,Xnorm) > tol)
  warning('CARE_NWT_FAC: no convergence in %d iterations\n', maxsteps);
end
