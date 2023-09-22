function [L, D, eqn, opts, oper] = ...
    expF(eqn, opts, oper, h, IQL, IQD, L_zero, D_zero, t_zero)
% Solve the affine problem arising from the split DRE.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin < 9
    t_zero = 0;
end

[out, eqn, opts, oper] = ...
    mess_exp_action(eqn, opts, oper, h, L_zero, t_zero);

L = [out.Z, IQL];
D = blkdiag(D_zero, IQD);

[L, D] = mess_column_compression(L, 'N', D, opts.splitting.trunc_tol, ...
                                 opts.splitting.trunc_info);
end
