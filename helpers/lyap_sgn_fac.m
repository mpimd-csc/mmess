function [R,iter] = lyap_sgn_fac(A,C,E)
%LYAP_SGN_FAC
%
% Solve the stable Lyapunov equation
%
% (*) A' X E + E' X A + C' C = 0
%                          X = R' R
%
% for a full-rank factor R of X via the sign function iteration.
%
% Input:
%    A - a square, n x n - matrix.
%    C - a p x n - matrix
%    E - a square, n x n - matrix.
%
% Output:
%    R - numerical full rank factor of X = R'*R.
%    iter (optionally) - number of sign iterations required.
%
% REFERENCE:
%
%   P. BENNER, E.S. QUINTANA-ORTI.
%   Model reduction based on spectral projection methods.
%   In: P. Benner, V.L. Mehrmann, D.C. Sorensen (eds.),
%   "Dimension Redution of Large-Scale Systems", vol. 45 of
%   Lecture Notes in Computational Science and Engineering, pp. 5-48,
%   Springer-Verlag, Berlin/Heidelberg, 2005.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


n = size(A,1);
R    = C;
if nargin < 3 || isempty(E)
    desc = 0;
    E    = eye(n);
    Enrm = 1;
else
    desc = 1;
    Enrm = norm(E,'fro');
end
Err     = norm(A + E,'fro');

% the following are parameters that could be used as input arguments in
% a pro version
maxstep = 50;
onemore = 0;
tol  = sqrt(n*eps)*Enrm;
rtol = 1e-8;

% further variables for convergence check
iter = 0;
convergence = Err <= tol;
In = eye(size(A,1));

while (iter < maxstep) && ((not(convergence)) || (convergence && (onemore < 2)))
    [AL,AU,p]=lu(A,'vector');
    p(p)=1:length(p);
    AL = AL( p , : );

    Y = AL\In;
    Y = AU\Y;
    if desc
        YE = Y*E;
        Y  = E*YE;
    end
    if Err > 0.1
        d = sqrt(norm(A,'fro')/norm(Y,'fro'));
    else
        d = 1;
    end
    A = (A/d + d*Y)/2;
    if desc
        R = [R; d*R*YE]/sqrt(2*d);
    else
        R = [R; d*R*Y]/sqrt(2*d);
    end

    [~,R,p]=qr(full(R),0);
    r = sum(abs(diag(R))>rtol*abs(R(1,1)));
    rc = size(R,2);
    q  = zeros(1,rc);
    for j=1:rc,  q(j) = find(p==j);  end
    R = R(1:r,q);

    Err  = norm(A + E,'fro');
    iter = iter + 1;
    convergence = Err <= tol;
    if convergence,  onemore = onemore + 1;  end
end
R = (R/E)/sqrt(2);
if (iter == maxstep) && (Err > tol)
     warning('LYAP_SGN_FAC: No convergence in %d iterations.\n', maxstep)
end
