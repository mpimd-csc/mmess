function X = lyap2solve(A,B)
% Solve Lyapunov equation AX+XA^T+B=0 via Zhou and Sorensen 2-solve
% method
%
% Usage:     X = lyap2solve(A,B)
%
% Input:
%  A         Matrix from Lyapunov equation
%  B         Matrix from Lyapunov equation
%
% Output:
%  X         Solution of Lyapunov equation

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


m = size(A,1);
n = size(B,2);

[Q,R]=schur(A);
idx = m:-1:1;

Q2=Q(:,idx);
R2=R(idx,idx)';

B = Q'*B*Q2;

Rsq = R*R;
I=speye(m);
j=1;

X=zeros(size(A));

while (j < n+1)

  if j==n ||...
          (j<n  && ...
           abs(R2(j+1,j)) < 10 * eps * max(abs(R2(j,j)), abs(R2(j+1,j+1))))

      if (j>1)
          b = -B(:,j) - X(:,1:j-1)*R2(1:j-1,j);
      else
          b = -B(:,j) ;
      end
      X(:,j) = (R+R2(j,j)*I)\b;
      j = j +1;
  else

      r11 = R2(j,j);  r12 = R2(j,j+1);
      r21 = R2(j+1,j); r22 = R2(j+1,j+1);

      if (j>1)
          b = -B(:,j:j+1) - X(:,1:j-1)*R2(1:j-1,j:j+1);
      else
          b = -B(:,j:j+1);
      end
      b = [R*b(:,1)+r22*b(:,1)-r21*b(:,2), R*b(:,2)+r11*b(:,2)-r12*b(:,1)];

      X(:,j:j+1) = ( Rsq+(r11+r22)*R + (r11*r22-r12*r21)*I)\b;

      j = j + 2;
  end
end

X = Q*X*Q2';


