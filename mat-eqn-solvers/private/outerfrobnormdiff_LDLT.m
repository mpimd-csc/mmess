function [s, s_ref] = outerfrobnormdiff_LDLT(L1, D1, L2, D2, rel)
% Compute the Frobenius norm of L1 D1 L1' - L2 D2 L2'
% for symmetric D1, D2

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

[L, D] = mess_column_compression([L1, L2], 'N', blkdiag(D1, -D2), ...
                                 eps, 0);
% Exploit the invariance under cyclic permutations property of the trace
% to compute norm(L*D*L','fro')
s = sqrt(trace((L' * L * D)^2));
if rel
    s_ref = sqrt(trace((L1' * L1 * D1)^2));
    s = s / s_ref;
else
    s_ref = [];
end

% Expensive version
%     s = norm(L1*D1*L1' - L2*D2*L2', 'fro');

end
