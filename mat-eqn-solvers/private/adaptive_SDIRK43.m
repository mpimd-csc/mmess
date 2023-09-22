function [out, eqn, opts, oper] = adaptive_SDIRK43(eqn, opts, oper, h, L, t_zero)
% Solve the system E(t)' \dot{z}(t) = A(t)' z(t) with z(0) = L over
% the interval [t0, t0 + h]. If t0 is omitted it is assumed to be 0.
% If eqn.type == 'N', instead solve E(t) \dot{z}(t) = A(t) z(t) .
% The method coefficients are taken from [1, p.100]. The method was chosen
% because it is of order 4 and therefore reasonably accurate in most
% situations. The SDIRK property further implies that its stability
% properties are much better than comparable explicit methods, while being
% much cheaper to evaluate than a fully implicit scheme. The embedded
% 3rd-order method provides the error estimate that allows for adaptive
% time-stepping.
%
% [1] Hairer, E. and Wanner, G., Solving Ordinary Differential Equations
% II: Stiff and Differential-Algebraic Problems, 2nd Ed., Springer, 2002.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

[n, p] = size(L);
normL = norm(L);
errest_order = 3;
tend = t_zero + h;
TOL = opts.exp_action.tol;

c = [1 / 4, 3 / 4, 11 / 20, 1 / 2, 1];

aij = [1 / 4,        0,        0,       0,    0
       1 / 2,       1 / 4,       0,       0,    0
       17 / 50,     -1 / 25,     1 / 4,      0,    0
       371 / 1360, -137 / 2720,  15 / 544,   1 / 4,   0
       25 / 24,    -49 / 48,   125 / 16,  -85 / 12, 1 / 4];

b = [25 / 24, -49 / 48, 125 / 16, -85 / 12, 1 / 4];

b_err = [-3 / 16, -27 / 32, 25 / 32, 0, 1 / 4];

onep = ones(1, p);

B = zeros(5 * p, 2 * p);

for k = 1:5

    B((k - 1) * p + 1:k * p, 1:p) = diag(b(k) * onep);
    B((k - 1) * p + 1:k * p, p + 1:2 * p) = diag(b_err(k) * onep);
end

C = zeros(5 * p, 5 * p);

for k = 1:5

    for l = 1:k - 1
        C((l - 1) * p + 1:l * p, (k - 1) * p + 1:k * p) = ...
            diag(aij(k, l) * onep);
    end
    C((k - 1) * p + 1:k * p, 1:p) =  diag(b(k) * onep);
end

t = t_zero;
hj = (tend - t_zero) / 100; % Initial step size, heuristic

if hj < eps
    out.Z = L;
    return
end

u = L;

% without eps the break condition can miss the last
% iteration step which leaves the field 'out' unset
tend_eff = tend + eps * tend;

while (t(end) + hj) < tend_eff

    K = zeros(n, 5 * p);

    for l = 1:5
        [eqn, opts, oper] = ...
            opts.splitting.eval_matrix_functions(eqn, opts, oper, ...
                                                 t(end) + c(l) * hj);
        RHS = oper.mul_A(eqn, opts, eqn.type, ...
                         u + hj * K(:, 1:(l - 1) * p) * ...
                         C(1:(l - 1) * p, (l - 1) * p + 1:l * p), 'N');

        % Want to do (M - aij*h*A)\RHS, but only have support for
        % (A + pE)\RHS, so scale by -1/(aih*h)
        coeff = -1.0 / (aij(l, l) * hj);

        K(:, (l - 1) * p + 1:l * p) = oper.sol_ApE(eqn, opts, eqn.type, ...
                                                   coeff, eqn.type, ...
                                                   coeff * RHS, 'N');
    end

    % 4th-order approximation:
    %     v2 = u + hj*(   25/24*k1 - 49/48*k2 +
    %                  + 125/16*k3 - 85/12*k4 + 1/4*k5);
    v = u + hj * K * B(:, 1:p);

    % 3rd-order approximation:
    %     v = u + h*(   59/48*k1 - 17/96*k2 +
    %                + 225/32*k3 - 85/12*k4);
    % Difference between 4th- and 3rd-order approximations,
    % evaluated in a better way:
    %     errest = hj*norm(  -3/16*k1 -27/32*k2 +
    %                      + 25/32*k3 + 0*k4 + 1/4*k5);
    errest = hj / (1 + normL) * norm(K * B(:, p + 1:2 * p));

    if errest > TOL % redo step
        hj = (0.9 * TOL / errest)^(1.0 / (errest_order)) * hj;
    else
        t(end + 1) = t(end) + hj; %#ok<AGROW>

        % New step size
        hj = (0.95 * TOL / errest)^(1.0 / errest_order) * hj;

        u = v;

        if abs(t(end) - tend) < eps % at final time
            out.Z = v;
            break
        end

        if t(end) + hj > tend % Ensure ending up precisely at tend
            hj = tend - t(end);
        end
    end
end
if not(exist('out', 'var'))
    out.Z = v;
    warning('something went wrong with the loop break');
end
end
