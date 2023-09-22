function [eqn, opts, oper] = KSM_mgs(eqn, opts, oper)
% KSM_MGS implements the (block) modified Gram-Schmidt
% orthogonalization method for the KSM framework
%
% Input and output:
%
%   eqn       structure containing equation data
%
%   opts      structure containing parameters for the algorithm
%
%   oper      contains function handles with operations for A and E
%
%
%   Most important here:
%   opts.KSM.compute_struct
%             structure that contain all the useful information computed in
%             the previous iterations with members:
%
%      opts.KSM.compute_struct.V     basis of the current space
%
%      opts.KSM.compute_struct.H     matrix collecting the
%                                          coefficients stemming from the
%                                          Arnoldi/Lanczos process
%
%      opts.KSM.compute_struct.T     projection of A onto the current
%                                          space, namely T=V'*A*V
%
%      opts.KSM.compute_struct.L,    matrices needed to recover the
%                                          columns of T
%      opts.KSM.compute_struct.rho   from the columns of H
%
%      opts.KSM.compute_struct.beta  p-by-p matrix such that
%                                          B=V_1*compute_struct.beta
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

p = size(eqn.W, 2);

VV = opts.KSM.compute_struct.V;
Vnew = opts.KSM.compute_struct.Vnew;
it = opts.KSM.compute_struct.it;
% perform a full orthogonalization for a better stability
k_max = it;

% current number of iteration
if strcmp('EK', opts.KSM.space)
    Hnew = zeros(2 * p * (it + 1), 2 * p);

    for l = 1:2
        k_min = max(1, it - k_max);
        for kk = k_min:it
            k1 = (kk - 1) * 2 * p + 1;
            k2 = kk * 2 * p;
            coef = VV(:, k1:k2)' * Vnew;
            Hnew(k1:k2, :) = Hnew(k1:k2, :) + coef;
            Vnew = Vnew - VV(:, k1:k2) * coef;
        end
    end

    % Normalization
    if it <= opts.KSM.maxiter
        [Vnew, Hnew(2 * p * it + 1:2 * p * (it + 1), :)] = qr(Vnew, 0);
    end

    % Compute also the projection of B if we are solving a care
    if strcmp(opts.KSM.type_eqn, 'CARE')
        switch eqn.type
            case 'N'
                Bnew = Vnew' * eqn.ssC';
            case 'T'
                Bnew = Vnew' * eqn.ssB;
        end
        opts.KSM.compute_struct.Bm = ...
            [opts.KSM.compute_struct.Bm; Bnew];
    end
    opts.KSM.compute_struct.V = [VV, Vnew];
    opts.KSM.compute_struct.H = Hnew;

elseif strcmp('RK', opts.KSM.space)
    H = opts.KSM.compute_struct.H;
    Hnew = zeros(p * (it + 1), p);

    for l = 1:2
        k_min = max(1, it - k_max);
        for kk = k_min:it
            k1 = (kk - 1) * p + 1;
            k2 = kk * p;
            coef = VV(:, k1:k2)' * Vnew;
            Hnew(k1:k2, :) = Hnew(k1:k2, :) + coef;
            Vnew = Vnew - VV(:, k1:k2) * coef;
        end
    end

    % Normalization
    if it <= opts.KSM.maxiter
        [Vnew, Hnew(p * it + 1:p * (it + 1), :)] = qr(Vnew, 0);
    end

    % Compute also the projection of B if we are solving a care
    if strcmp(opts.KSM.type_eqn, 'CARE')
        switch eqn.type
            case 'N'
                Bnew = Vnew' * eqn.ssC';
            case 'T'
                Bnew = Vnew' * eqn.ssB;
        end
        opts.KSM.compute_struct.Bm = ...
            [opts.KSM.compute_struct.Bm; Bnew];
    end
    opts.KSM.compute_struct.V = [VV, Vnew];

    if not(isempty(H))
        opts.KSM.compute_struct.H = [[H; zeros(p, p * (it - 1))], Hnew];
    else
        opts.KSM.compute_struct.H = Hnew;
    end
end

end
