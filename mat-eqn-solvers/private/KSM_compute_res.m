function [eqn, opts, oper] = KSM_compute_res(eqn, opts, oper)
% function [eqn, opts, oper] = KSM_compute_res(eqn,opts,oper)
% Function that compute the absolute residual norm
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
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

it = opts.KSM.compute_struct.it;

if not(isfield(opts.KSM.compute_struct, 'res'))
    opts.KSM.compute_struct.res = [];
end

if strcmp('EK', opts.KSM.space)
    p = size(opts.KSM.compute_struct.beta, 2) / 2;
    Y = opts.KSM.compute_struct.Y;
    T = opts.KSM.compute_struct.T;

    cc = T(2 * p * it + 1:2 * p * (it + 1), end - 2 * p + 1:end);
    % Compute the residual norm
    if eqn.haveE
        W = [opts.KSM.compute_struct.V(:, 1:2 * p * it) * ...
             (Y(:, 2 * p * (it - 1) + 1:2 * p * it) * cc'), ...
             opts.KSM.compute_struct.V(:, 2 * p * it + 1: ...
                                       2 * p * (it + 1))];
        switch eqn.type
            case 'N'
                D = oper.ss_to_dss(eqn, opts, 'L', 'N', W, 'N');
            case 'T'
                D = oper.ss_to_dss(eqn, opts, 'U', 'T', W, 'N');
        end

        [~, N] = qr(full(D), 0);
        N = triu(N);
        K = zeros(4 * p);
        K(1:2 * p, 2 * p + 1:end) = eye(2 * p);
        K(2 * p + 1:end, 1:2 * p) = eye(2 * p);
        res = norm(N * K * N', 'fro');
    else
        res = sqrt(2) * ...
            norm(cc * Y(2 * p * (it - 1) + 1:2 * p * it, :), 'fro');
    end
    opts.KSM.compute_struct.res = ...
        [opts.KSM.compute_struct.res, res];

elseif strcmp('RK', opts.KSM.space)
    p = size(opts.KSM.compute_struct.beta, 2);
    Y = opts.KSM.compute_struct.Y;
    T = opts.KSM.compute_struct.T;
    H = opts.KSM.compute_struct.H;
    V = opts.KSM.compute_struct.V;

    W = [V(:, 1:p * it) * ...
         (Y * (H(1:p * it, 1:p * it)' \ ...
               H(p * it + 1:p * (it + 1), 1:p * it)')), ...
         V(:, p * it + 1:p * (it + 1)) * ...
         opts.KSM.compute_struct.shifts(end - 1) - ...
         opts.KSM.compute_struct.AVnew + ...
         V(:, 1:p * it) * T(1:p * it, p * it + 1:p * (it + 1))];

    if eqn.haveE
        switch eqn.type
            case 'N'
                D = oper.ss_to_dss(eqn, opts, 'L', 'N', W, 'N');
            case 'T'
                D = oper.ss_to_dss(eqn, opts, 'U', 'T', W, 'N');
        end
    else
        D = W;
    end
    [~, N] = qr(full(D), 0);
    N = triu(N);
    K = zeros(2 * p);
    K(1:p, p + 1:end) = eye(p);
    K(p + 1:end, 1:p) = eye(p);
    res = norm(N * K * N', 'fro');
    opts.KSM.compute_struct.res = ...
        [opts.KSM.compute_struct.res, res];

end

end
