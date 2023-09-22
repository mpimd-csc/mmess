function eqn = LU_E(eqn)
% LU_E computes and caches LU decomposition of eqn.E_

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if isfield(eqn, 'haveE') && eqn.haveE
    if isfield(eqn, 'EL') && isfield(eqn, 'EU') && ...
            isfield(eqn, 'Ep') && isfield(eqn, 'Eq') && ...
            isfield(eqn, 'ER')
        if isfield(eqn, 'Ecount')
            eqn.Ecount = eqn.Ecount + 1;
        else
            eqn.Ecount = 2;
        end
    else
        if issymmetric(eqn.E_)
            % try to compute the Cholesky factorization if E is symmetric
            [eqn.EL, p, eqn.Ep] = chol(eqn.E_, 'lower', 'vector');
        else
            p = 1;
        end
        % if E is also positive definite, p==0
        if p == 0
            eqn.EU = eqn.EL';
            eqn.Eq = eqn.Ep;
            eqn.ER = speye(size(eqn.E_));
            % reorder eqn.ER already once so that when we reapply
            % the reordering in
            % ss_to_dss_state_space_transformed_default we get the identity
            % eqn.ER = eqn.ER(eqn.Ep,:);
        else
            [eqn.EL, eqn.EU, eqn.Ep, eqn.Eq, eqn.ER] = lu(eqn.E_, 'vector');
        end
        eqn.iEq(eqn.Eq) = 1:size(eqn.E_, 1);
        eqn.Ecount = 1;
    end
end
