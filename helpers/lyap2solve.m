function X = lyap2solve(A, B)
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
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

m = size(A, 1);
n = size(B, 2);

[Q, R] = schur(A);
idx = m:-1:1;

Q2 = Q(:, idx);
R2 = R(idx, idx)';

B = Q' * B * Q2;

Rsq = R * R;
Id = speye(m);
k = 1;

X = zeros(size(A));

tol = 10 * eps;

while k < (n + 1)

    if k == n || ...
       (k < n  && ...
        abs(R2(k + 1, k)) < tol * max(abs(R2(k, k)), ...
                                      abs(R2(k + 1, k + 1))))

        if k > 1
            b = -B(:, k) - X(:, 1:k - 1) * R2(1:k - 1, k);
        else
            b = -B(:, k);
        end

        X(:, k) = (R + R2(k, k) * Id) \ b;
        k = k + 1;
    else

        r11 = R2(k, k);
        r12 = R2(k, k + 1);
        r21 = R2(k + 1, k);
        r22 = R2(k + 1, k + 1);

        if k > 1
            b = -B(:, k:k + 1) - X(:, 1:k - 1) * R2(1:k - 1, k:k + 1);
        else
            b = -B(:, k:k + 1);
        end

        b = [R * b(:, 1) + r22 * b(:, 1) - r21 * b(:, 2), ...
             R * b(:, 2) + r11 * b(:, 2) - r12 * b(:, 1)];

        X(:, k:k + 1) = (Rsq + (r11 + r22) * R + ...
                         (r11 * r22 - r12 * r21) * Id) \ b;

        k = k + 2;
    end
end

X = mess_symmetrize(Q * X * Q2');
