function [TL, TR] = square_root_method_SO(M, max_ord, tol, inform, U, V)
% Square root method for the computation of the transformation matrices to
% balance and reduce second order systems
%
% Call
%
%  [TL,TR] = square_root_method_SO(M, max_ord, tol, inform, U,V)
%
% Inputs:
%  M                the second order systems mass matrix
%
%  max_ord          maximum allowed reduced order
%
%  tol              truncation tolerance
%                   (drops all singular values smaller than
%                    tol*largest singular value)
%
%  ZB, ZC           the (tall and skinny) Gramian factors
%
% Outputs:
%  TL,TR            left and right truncation matrices
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

[U0, S0, V0] = svd(V' * M * U, 0);
opts = struct;

s0 = diag(S0);
ks = length(s0);
k0 = 1;
while (k0 <= ks) && (s0(k0) > s0(1) * tol)
    k0 = k0 + 1;
end

r = min([max_ord k0]);
if inform > 0
    mess_fprintf(opts, ['reduced system order: %d ', ...
                        '(max possible/allowed: %d/%d)\n\n'], r, ks, max_ord);
end

sigma_r = diag(S0(1:r, 1:r));

VB = U * V0(:, 1:r);
VC = V * U0(:, 1:r);

TR = VB * diag(ones(r, 1) ./ sqrt(sigma_r));
TL = VC * diag(ones(r, 1) ./ sqrt(sigma_r));
end
