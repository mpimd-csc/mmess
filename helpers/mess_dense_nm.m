function X = mess_dense_nm(opts, A, B, C, E, X0, Q, R)
% naive Newton Kleinman iteration for the ARE
%
%     A'*X*E + E'*X*A + C'*Q*C - E'*X*B*R\B'*X*E = 0
%
% Inputs:
%
% opts              options structure needed for the logger
% A, B, C, E, Q, R  Coefficients in the above equation
% X0                initial guess for the solution
%
% Outputs:
% X                 Solution
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%%
% check for  and select available Lyapunov solver
if exist('lyap', 'file')
    meth = 'lyap';

elseif exist('lyap2solve', 'file')
    % Again shipped with our code so actually unlikely to be unavailable.
    meth = 'lyap2solve';

elseif exist('lyap_sgn_fac', 'file')
    % This should always be available since we ship it.
    % The following are only here in case someone deletes the sign solver
    meth = 'lyap_sgn_fac';

elseif exist('lyapchol', 'file')
    % In case lyap is not available it is rather unlikely that this will be
    % there. Still someone might have their own implementation...
    meth = 'lyapchol';

else
    mess_err(opts, 'missing_solver', ...
             'mess_dense_nm was unable to find Lyapunov solver');

end

%%
% prepare data for the constant term
if (nargin > 6) && not(isempty(Q))
    if not(issymmetric(Q))
        mess_err(opts, 'check_data', 'Q must be symmetric');
    end

    % Q must be symmetric pos. semidef. for lyapchol or lyap_sgn_fac
    if strcmp(meth, 'lyapchol') || strcmp(meth, 'lyap_sgn_fac')
        [U, S_Q] = eig(Q);
        if any(diag(S_Q) < 0)
            meth = 'lyap2solve';
        end
    end

    G = C' * Q * C;
    switch meth
        case 'lyap'
            G = (G + G') / 2; % make sure it's symmetric for e.g. lyap

        case {'lyapchol', 'lyap_sgn_fac'}
            C = sqrt(S_Q) * U' * C; % C'*C = G

        case 'lyap2solve'
            if (nargin < 4) || isempty(E)
                GE = C' * Q * C;

            else
                CE = C / E;
                GE = CE' * Q * CE;

            end
    end

else

    G = C' * C;
    switch meth
        case 'lyap'
            G = (G + G') / 2; % make sure it's symmetric for e.g. lyap

        case 'lyap2solve'
            if (nargin < 5) || isempty(E)
                GE = C' * C;

            else
                CE = C / E;
                GE = CE' * CE;

            end
    end

end

res0 = norm(G);

%%
% prepare center matrix in the quadratic term
if (nargin < 5) || isempty(E)
    E = eye(size(A, 1));
end
if (nargin > 7) && not(isempty(R))
    if not(issymmetric(R))
        mess_err(opts, 'check_data', 'R must be symmetric');
    end
    F = B * (R \ B');
else
    F = B * B';
    R = eye(size(B, 2));
end

%% Main Newton loop
tol = 1e-14;
maxiter = 50;

for k = 1:maxiter
    if k > 1 || ((nargin > 5) && not(isempty(X0)))
        K = R \ (B' * X0 * E);
    else
        K = zeros(size(B'));
        X0 = zeros(size(A));
    end
    switch meth
        case 'lyap'
            RHS = G + K' * R * K;
            RHS = .5 * (RHS + RHS');
            X = lyap(A' - K' * B', RHS, [], E');
        case 'lyapchol'
            XC = lyapchol(A' - K' * B', [C', K' * sqrtm(R)], E');
            X = XC' * XC;
        case 'lyap_sgn_fac'
            XC = lyap_sgn_fac(A - B * K, [C; K], E);
            X = XC' * XC;
        case 'lyap2solve'
            KE = K / E;
            X = lyap2solve(((A - B * K) / E)', GE + KE' * R * KE);
    end
    XE = X * E;
    res = norm(A' * XE + XE' * A - XE' * F * XE + G);
    rc = norm(X - X0) / norm(X);

    if (rc < tol) || (res < tol * res0)
        break
    else
        X0 = X;
    end
end

%% Let us warn the user if we stopped before converging
if k == maxiter && not(rc < tol) && not(res < tol * res0)
    mess_warn(opts, 'denseNM_convergence', ...
              ['dense Newton method stopped by maximum iteration number.', ...
               ' Results may be inaccurate!']);
end
