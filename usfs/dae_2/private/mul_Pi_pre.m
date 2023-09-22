function [eqn, opts, oper] = mul_Pi_pre(eqn, opts, oper)
% MUL_Pi multiplies with the hidden manifold projection matrix or it
% transpose. Note that the multiplication is actually implemented as the
% solution of a saddle point system.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isfield(eqn, 'P_'))
    eqn.P_ = eqn.A_;
    eqn.P_(1:eqn.manifold_dim, 1:eqn.manifold_dim) = eqn.E_(1:eqn.manifold_dim, 1:eqn.manifold_dim);
    eqn.Pcount = 1;
else
    eqn.Pcount = eqn.Pcount + 1;
end
