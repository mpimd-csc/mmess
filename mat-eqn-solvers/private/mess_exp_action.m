function [out, eqn, opts, oper] = ...
    mess_exp_action(eqn, opts, oper, h, L, t_zero)
% Computes the matrix exponential action expm(h*(E\A))*L where L is a
% skinny block matrix.
%
% Computes (expm(h*(E'\A'))*L if eqn.type == 'T'
% Computes (expm(h*(E\A))*L if eqn.type == 'N'
%
% If t0 is given as input and eqn.LTV == true, instead solve the LTV system
%   E(t)'\dot{x}(t) = A(t)' x(t) (eqn.type == 'T')
%   E(t)\dot{x}(t) = A(t) x(t)   (eqn.type == 'N')
% over the interval [t0, t0 + h].
%
% NOTE: Only the Krylov and adaptive SDIRK43 methods are implemented for
% the autonomous case so far! Only adaptive SDIRK43 for the LTV case,
% through 'LTV'.
%
%
% Input and output:
%   eqn       structure containing equation data
%
%   opts      structure containing parameters for the algorithm
%
%   oper      contains function handles with operations for A and E
%
% Input:
%
%   h         contains step size h
%
%   L         contains block matrix L
%
%   t_zero    contains starting time t0 in the LTV case
%
% Input fields in opts.exp_action:
%
%   opts.exp_action.method
%       possible values:  'Krylov', 'adaptiveSDIRK', 'LTV'
%       use this method to compute matrix exponential actions on block
%       matrices
%       (optional, default 'Krylov')
%
%   opts.exp_action.tol
%       possible values: scalar > 0
%       tolerance for matrix exponential actions (means different things
%       for Krylov and adaptiveSDIRK)
%       (optional, default 1e-4)
%
%   opts.exp_action.Krylov.kabsmax
%       possible values: positive integer
%       the absolute maximum of Krylov iterations to run (due to memory
%       limitations)
%       (optional, defaults to a number corresponding to about 4GB memory)
%
%
% Output:
%  out              structure containing the following:
%
%  out.Z            the matrix exponential action
%
%  out.converged    flag for the iterative methods, which is 1 if the
%                   method converged and 0 otherwise
%
%  out.errest       the final error estimate (residual)

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for matrix exponential actions control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isfield(opts, 'exp_action')) || ...
        not(isstruct(opts.exp_action)) || ...
        not(isfield(opts.exp_action, 'method'))
    mess_warn(opts, 'control_data', ...
              ['matrix exponential actions control ' ...
               'structure opts.exp_action ' ...
               'missing.  Using default Krylov method.']);
    opts.exp_action.method = 'Krylov';
end

if not(isfield(opts.exp_action, 'tol'))
    opts.exp_action.tol = 1e-4;
end

if strcmp(opts.exp_action.method, 'Krylov')
    if eqn.LTV
        mess_err(opts, 'control_data', ...
                 ['LTV problem specified, but opts.exp_action.method=', ...
                  '''Krylov''. Use option ''LTV'' instead.']);
    end

    if not(isfield(opts.exp_action, 'Krylov')) || ...
        not(isstruct(opts.exp_action.Krylov))
        opts.exp_action.Krylov = {};
    end

    if not(isfield(opts.exp_action.Krylov, 'kabsmax'))
        [n, p] = size(L);

        mmem = 8 * 4 * 1024^3; % 4 GB
        kabsmax = mmem / n / p;
    else
        kabsmax = opts.exp_action.Krylov.kabsmax;
    end

elseif strcmp(opts.exp_action.method, 'adaptiveSDIRK')
    % No extra parameters for this method. Tolerance
    % is set globally in opts.exp_action.tol
    if eqn.LTV
        mess_err(opts, 'control_data', ...
                 ['LTV problem specified, but opts.exp_action.method=', ...
                  '''adaptiveSDIRK''. Use option ''LTV'' instead.']);
    end

    if not(eqn.LTV)
        % Assume time interval [0, h] for backwards compatibility
        if nargin < 6
            t_zero = 0;
        end
    end

elseif strcmp(opts.exp_action.method, 'LTV')
    % Like the previous method, this one also needs no extra parameters
    if not(eqn.LTV)
        mess_err(opts, 'control_data', ...
                 'opts.exp_action.method = ''LTV'' but eqn.LTV = false.');
    end

    if nargin < 6  % No t0 given
        mess_err(opts, 'control_data', ...
                 ['LTV problem specified, but no initial time set for ' ...
                  'call to mess_exp_action']);
    end

else
    mess_err(opts, 'control_data', ...
             ['Chosen method for matrix exponential ' ...
              'actions: ', opts.exp_action.method, ' is not supported']);

end

switch opts.exp_action.method

    case 'Krylov'
        [~, p] = size(L);
        normL = norm(L);
        tol = opts.exp_action.tol;

        Afun = @(x) -h * oper.sol_E(eqn, ...
                                    opts, ...
                                    eqn.type, ...
                                    oper.mul_A(eqn, ...
                                               opts, ...
                                               eqn.type, ...
                                               x, ...
                                               'N'), ...
                                    'N');

        kmin = 1;
        kmax = 2; % Do two initial blocks before estimating the error

        [Vk, R] = qr(L, 0); % Here, Vk = U1
        Hk = [];

        converged = false;
        while not(converged)

            for k = kmin:kmax
                Uk = Vk(:, (k - 1) * p + 1:k * p);

                Wk = Afun(Uk);

                for i = 1:k % Orthogonalize
                    Ui = Vk(:, (i - 1) * p + 1:i * p);
                    Hik = Ui' * Wk;
                    Wk = Wk - Ui * Hik;
                    Hk((i - 1) * p + 1:i * p, (k - 1) * p + 1:k * p) = Hik;
                end

                [Ukp1, Hkp1] = qr(Wk, 0);

                if k == kmax
                    % Compute error estimate
                    E1 = eye(k * p, k * p);
                    E1 = E1(:, 1:p);
                    Hkt = [-Hk, E1; zeros(p, (k + 1) * p)];

                    eHt = expm(Hkt);
                    phiHkE1 = eHt(1:k * p, k * p + 1:(k + 1) * p);

                    errest = normL * norm(Hkp1) * ...
                             norm(phiHkE1((k - 1) * p + 1:k * p, :));

                    if errest < tol
                        out.converged = true;
                        out.errest = errest;
                        % The below essentially does eHk = expm(-Hk) and
                        % Z =  Vk*expm(-Hk)*(Vk'*L), but in a more
                        % efficient way
                        eHk = eHt(1:k * p, 1:k * p);
                        out.Z = Vk * (eHk(:, 1:p) * R);
                        return
                    end
                end

                Vk(:, k * p + 1:(k + 1) * p) = Ukp1;
                Hk(k * p + 1:(k + 1) * p, (k - 1) * p + 1:k * p) = Hkp1;

            end

            kmin = kmax + 1;
            kmax = kmax + 3; % Add two blocks in every step

            if kmax > kabsmax
                out.Z = NaN; % TODO: improve error handling, see #312
                out.converged = false;
                out.errest = Inf;
                mess_warn(opts, 'exp_action', ...
                          ['Krylov method for matrix exponential ' ...
                           'action did NOT converge!']);
                return
            end

        end

    case 'adaptiveSDIRK'
        [out, eqn, opts, oper] = ...
            adaptive_SDIRK43(eqn, opts, oper, h, L, t_zero);

    case 'LTV'
        % Temporarily change the matrix updating function to evaluate at
        % t0 + h - s instead of at s
        opts.splitting.eval_matrix_functions_temp = ...
            opts.splitting.eval_matrix_functions;

        opts.splitting.eval_matrix_functions = @(eqn, opts, oper, s) ...
            opts.splitting.eval_matrix_functions_temp(eqn, opts, oper, ...
                                                      t_zero + h - s);

        [out, eqn, opts, oper] = adaptive_SDIRK43(eqn, opts, oper, ...
                                                  h, L, 0);

        % Restore the matrix updating function
        opts.splitting.eval_matrix_functions = ...
            opts.splitting.eval_matrix_functions_temp;
end
