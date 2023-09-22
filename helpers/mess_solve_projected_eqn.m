function [eqn, opts, oper] = ...
    mess_solve_projected_eqn(eqn, opts, oper, framework, type)
% function [out, eqn, opts, oper] =
%       mess_solve_projected_eqn(eqn, opts, oper, framework, type)
%
% Function that solves small-scale projected Lyapunov and Riccati
% equations, either projected to the solution span (Galerkin projection
% acceleration (GPA), e.g. in LRNM) or the Krylov basis in KSM.
%
% Input and output:
%
%  type                     possible values: 'LE','CARE'
%                           determines whether a Lyapunov ('LE') or a
%                           Riccati ('CARE') equation should be
%                           projected
%
%  eqn                      structure with data for A, E, W
%                           in the equation determined by type
%
%  oper                     structure contains function handles for
%                           operations with A, E
%
%  framework                possible values: 'KSM','GPA'
%                           KSM: the projected equation involved in the
%                                Krylov method is solved
%                           GPA: the Galerkin acceleration scheme is
%                                applied
%
%  opts                     options structure that needs to contain one
%                           of the following members
%
%  opts.adi                 options structure for ADI method
%
%  opts.nm                  options structure for Newton method
%
%  opts.KSM                 options structure for Krylov method
%
%
%  xopts                    opts.adi, opts.nm or opts.KSM depending on
%                           'type' above. We expect a substructure called
%                           'projection' in xopts that has entries:
%
%     ortho                 possible  values: false, true
%                           implicit (false) or explicit (true)
%                           orthogonalization of Z in LRNM
%                           i.e., explicit orthogonalization via orth().
%                           (optional, default: true,
%                            ignored in framework KSM)
%
%     meth                  method for solving projected Lyapunov or
%                           Riccati equation. Depending on 'type' possible
%                           values are: 'lyapchol', 'lyap_sgn_fac', 'lyap',
%                           'lyapunov', 'lyap2solve', 'icare', 'care',
%                           'care_nwt_fac', 'mess_dense_nm'
%                           (optional, default: best available solver for the
%                            type depending on the presence of the control
%                            toolbox/package and version of Matlab/Octave used.)
%                           Remark: some solver are disallowed/excluded when
%                           opts.LDL_T is 'true'.
%
% Output:
%
% (KSM)
% opts.KSM.compute_struct.Y  solution of the projected equation
%
% (GPA)
% xopts.projection.Z
% xopts.projection.D       Updated solution factors after prolongation written
%                          into the correct substructure xopts represents
%
%
% uses operatorfunctions mul_A, mul_E, mul_ApE

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check Inputs
if not(isfield(opts, 'LDL_T')) || isempty(opts.LDL_T)
    opts.LDL_T = false;
end

switch framework
    case 'GPA'

        switch type

            case {'LE'}

                xopts = opts.adi;

            case {'CARE'}

                xopts = opts.nm;

            otherwise

                mess_err(opts, 'GP_type', ...
                         ['type has to be ''LE'' or ''CARE'' in Galerkin '...
                          'projection acceleration.']);
        end

    case 'KSM'

        xopts = opts.KSM;

    otherwise
        mess_err(opts, 'GP_framework', ...
                 ['unknown framework requested. Must be either', ...
                  ' ''GPA'' or ''KSM''.']);
end

% check for process required parameters in xopts
if not(isfield(xopts, 'projection')) || ...
        not(isfield(xopts.projection, 'ortho')) || ...
        not(islogical(xopts.projection.ortho)) % last one covers non-empty

    switch framework
        case 'GPA' % we want to orthogonalize the factor columns

            xopts.projection.ortho = true;

        case 'KSM' % orthogonalization is part of the basis generation

            xopts.projection.ortho = false;
    end
end

% if opts.projection.meth is not set, or the user did not specify an actual
% solver method but an equation let us set a default solver
if not(isfield(xopts.projection, 'meth')) || ...
        isempty(xopts.projection.meth) || ...
        isequal(xopts.projection.meth, 'lyapunov') || ...
        isequal(xopts.projection.meth, 'riccati')

    switch type
        case 'LE'
            xopts.projection.meth = find_recommended_CALE_solver(opts.LDL_T);
        case 'CARE'
            xopts.projection.meth = find_recommended_CARE_solver();
    end
end

% now we check if the requested routines are actually available and set a
% fallback if not.
if not(exist(xopts.projection.meth, 'file'))
    meth_non_ex = xopts.projection.meth;
    switch type
        case 'LE'
            xopts.projection.meth = find_recommended_CALE_solver(opts.LDL_T);
        case 'CARE'
            xopts.projection.meth = find_recommended_CARE_solver();
    end

    mess_warn(opts, 'missing_solver', ...
              ['mess_solve_projected_eqn was unable to find'...
               ' solver ''%s'', switched to ''%s'''], meth_non_ex, ...
              xopts.projection.meth);
end

% In the LDL_T case the core matrix may be indefinite and ZZ^T factored dense
% solvers may fail, or use complex data that we carefully avoid otherwise, so we
% disallow them
if opts.LDL_T && ...
        (isequal(xopts.projection.meth, 'care_nwt_fac') || ...
         isequal(xopts.projection.meth, 'lyap_sgn_fac') || ...
         isequal(xopts.projection.meth, 'lyapchol'))

    mess_err(opts, 'illegal_input', ...
             ['We do not allow ZZ^T factored dense solver in LDL_T mode,'...
              'as ''D'' may be indefinite.']);

end

%% Set small matrices

switch framework
    case 'KSM'
        % prepare the coefficients of the projected equation
        beta = opts.KSM.compute_struct.beta;
        it = opts.KSM.compute_struct.it;
        A = opts.KSM.compute_struct.T;
        if strcmp('EK', opts.KSM.space)
            p = size(beta, 2) / 2;
            A = A(1:2 * p * it, 1:2 * p * it);
            B = eye(size(A, 2), p) * beta(1:p, 1:p) * ...
                diag(sqrt(diag(eqn.T)));
        elseif strcmp('RK', opts.KSM.space)
            p = size(beta, 2);
            A = A(1:p * it, 1:p * it);
            B = eye(size(A, 2), p) * beta(1:p, 1:p) * ...
                diag(sqrt(diag(eqn.T)));
            S = 1.0; % actually eye(size(B, 2)); but scalar is cheaper
            U = 1.0; % actually eye(size(B, 2)); but scalar is cheaper

        end

        if strcmp(type, 'CARE')
            A = A';
            C = B';
            S = 1.0; % actually eye(size(C, 1)); but scalar is cheaper
            U = 1.0; % actually eye(size(C, 1)); but scalar is cheaper
            B = opts.KSM.compute_struct.Bm;
            if strcmp('EK', opts.KSM.space)
                B = B(1:2 * p * it, :);
            elseif strcmp('RK', opts.KSM.space)
                B = B(1:p * it, :);
            end
            R = eye(size(B, 2));
        end
        if exist('OCTAVE_VERSION', 'builtin')
            E = eye(size(A));
        else
            E = [];
        end

    case 'GPA'
        if not(xopts.projection.ortho)
            mess_warn(opts, 'accuracy', ...
                      ['Galerkin projection acceleration without ', ...
                       'orthogonalization of the projection basis ', ...
                       'is dangerous. Only use this when you know that ', ...
                       'orthogonalization happened outside already.']);
        end
        [A, B, C, E, Z, S, U, R, lyapunov] = ...
            GPA_projected_matrices(eqn, opts, oper, xopts, type);
end

%% Here comes the actual solution process for the projected matrix equation
factorize = true;
switch type
    case 'LE'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Choose solver for the small equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch xopts.projection.meth
            case 'lyapchol'
                if opts.LDL_T
                    B = B * U * sqrt(S);
                    if isempty(E)
                        XC = lyapchol(A, B);
                    else
                        XC = lyapchol(A, B, E);
                    end
                    [~, S, XC] = svd(XC, 'econ');
                    XC = XC';
                    xopts.projection.D = S.^2;
                else
                    if isempty(E)
                        XC = lyapchol(A, B);
                    else
                        XC = lyapchol(A, B, E);
                    end
                end
                factorize = false;

            case 'lyap_sgn_fac'
                if opts.LDL_T
                    B = B * U * sqrt(S);
                    XC = lyap_sgn_fac(A', B', E');
                    [~, S, XC] = svd(XC, 'econ');
                    XC = XC';
                    xopts.projection.D = S.^2;
                else
                    XC = lyap_sgn_fac(A', B', E');
                    XC = XC';
                end
                factorize = false;

            case {'lyap'}
                if opts.LDL_T
                    B = B * U * S * U' * B';
                    B = (B + B') / 2; % make sure it's symmetric for lyap
                    if isempty(E)
                        X = lyap(A, B);
                    else
                        X = lyap(A, B, [], E);
                    end
                else
                    if isempty(E)
                        X = lyap(A, B * B');
                    else
                        X = lyap(A, B * B', [], E);
                    end
                end

            case 'lyap2solve'
                if eqn.haveE && not(isempty(E))
                    EB = E \ B;
                    if opts.LDL_T
                        X = lyap2solve(E \ A, EB * U * S * U' * EB');
                    else
                        X = lyap2solve(E \ A, EB * EB');
                    end
                else
                    if opts.LDL_T
                        X = lyap2solve(A, B * U * S * U' * B');
                    else
                        X = lyap2solve(A, B * B');
                    end

                end
        end

    case 'CARE'
        switch xopts.projection.meth
            case {'care'}
                if opts.LDL_T
                    G = C' * S * C;
                else
                    G = C' * C;
                    S = [];
                end
                if isempty(E)
                    X = care(A, B, G, R);
                    res = norm(A' * X + X * A - ...
                               (X * B) * (R \ (B' * X)) +  G, 'fro');
                else
                    X = care(A, B, G, R, [], E);
                    res = norm(A' * X * E + E' * X * A - ...
                               E' * (X * B) * (R \ (B' * X)) * E +  G, 'fro');
                end
                if res / norm(G, 'fro') >= 1e-12
                    X = mess_dense_nm(opts, A, B, C, E, X, S, R);
                end
            case {'icare'}
                if opts.LDL_T
                    G = C' * S * C;
                else
                    G = C' * C;
                    S = [];
                end
                if isempty(E)
                    X = icare(A, B, G, R);
                    res = norm(A' * X + X * A - ...
                               (X * B) * (R \ (B' * X)) +  G, 'fro');
                else
                    X = icare(A, B, G, R, [], E);
                    res = norm(A' * X * E + E' * X * A - ...
                               E' * (X * B) * (R \ (B' * X)) * E +  G, 'fro');
                end
                if res / norm(G, 'fro') >= 1e-12
                    X = mess_dense_nm(opts, A, B, C, E, X, S, R);
                end
            case 'care_nwt_fac'
                if opts.LDL_T
                    C = sqrt(S) * U' * C;
                end
                if not(isempty(E))
                    XC = care_nwt_fac([], A / E, B, C / E, 1e-10, 100);
                    X = XC' * XC;
                    if norm(A' * X * E + E' * X * A - ...
                            (E' * (X * B)) * (B' * X) * E + ...
                            C' * S * C, 'fro') / ...
                       norm(C' * S * C, 'fro') >= 1e-12
                        X = mess_dense_nm(opts, A, B, C, [], X, S, R);
                    end
                else
                    XC = care_nwt_fac([], A, B, C, 1e-10, 100);
                    X = XC' * XC;
                    if norm(A' * X + X * A - ...
                            (X * B) * (B' * X) + ...
                            C' * S * C, 'fro') / ...
                       norm(C' * C, 'fro') >= 1e-12
                        X = mess_dense_nm(opts, A, B, C, [], X, S, R);
                    end
                end
                if opts.LDL_T
                    [~, S, XC] = svd(XC, 'econ');
                    XC = XC';
                    opts.projection.D = diag(S).^2;
                    XC = XC';
                else
                    clear XC;
                end
                factorize = false;

            case 'mess_dense_nm'
                if opts.LDL_T
                    X = mess_dense_nm(opts, A, B, C, E, [], S, R);
                else
                    X = mess_dense_nm(opts, A, B, C, E);
                end
        end

end
%% Post processing
% For KSM we actually want the full solution matrix, while for GPA need some
% factored form of it.
switch framework

    case 'KSM'
        if exist('XC', 'var')
            if opts.LDL_T
                X = XC * xopts.projection.D * XC';
            else
                X = XC * XC';
            end
        end
        opts.KSM.compute_struct.Y = X;

    case 'GPA'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % If the projected solution was not already computed in factored
        % form compute a symmetric factorization now and update the large
        % factor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if factorize
            if opts.LDL_T

                xopts.projection.D = X;
                xopts.projection.Z = Z;

            elseif exist('cholp', 'file')
                [XC, P, I] = cholp(X);
                XC = P * XC';

                if I && lyapunov
                    mess_warn(opts, 'proj_sol_semidef', ...
                              ['The solution of the projected ', ...
                               'equation was semidefinite.']);
                end
            else
                [~, S, V] = svd(X);
                s = diag(S);
                r = find(s > s(1) * eps);
                xopts.projection.XC = diag(sqrt(s(r))) * V(:, r)';
            end
        end

        if exist('XC', 'var')
            xopts.projection.Z = Z * XC;
        end

end

%% Finally we need to write the results to the correct substructure
switch framework
    case 'GPA'

        switch type

            case {'LE'}

                opts.adi.projection = xopts.projection;

            case {'CARE'}

                opts.nm.projection = xopts.projection;

            otherwise

                mess_err(opts, 'GP_type', ...
                         ['type has to be ''LE'' or ''CARE'' in Galerkin '...
                          'projection acceleration.']);
        end

    case 'KSM'

        opts.KSM.projection = xopts.projection;

end

end % of main function

%% Local helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the recommended solver for the continuous ALE
function out = find_recommended_CALE_solver(LDL_T)
if exist('lyap', 'file')
    out = 'lyap';
elseif LDL_T && exist('lyap2solve', 'file')
    out = 'lyap2solve';
elseif exist('lyap_sgn_fac', 'file')
    out = 'lyap_sgn_fac';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the recommended solver for continuous ARE
function out = find_recommended_CARE_solver()

if exist('icare', 'file')
    out = 'icare';
elseif exist('care', 'file')
    out = 'care';
elseif exist('mess_dense_nm', 'file')
    out = 'mess_dense_nm';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute projected coefficient matrices
function [A, B, C, E, Z, S, U, R, lyapunov] = ...
    GPA_projected_matrices(eqn, opts, oper, xopts, type)

if not(isfield(opts, 'rosenbrock'))
    opts.rosenbrock = [];
end
if isstruct(opts.rosenbrock) && isfield(opts.rosenbrock, 'tau')
    rosenbrock = 1;
    if opts.rosenbrock.stage == 1
        pc = -1 / (2 * opts.rosenbrock.tau);
    else % stage 2
        pc = -1 / (2 * opts.rosenbrock.tau * opts.rosenbrock.gamma);
    end
else
    rosenbrock = 0;
end
if not(isfield(opts, 'bdf'))
    opts.bdf = [];
end
if isstruct(opts.bdf) && isfield(opts.bdf, 'tau') && ...
        isfield(opts.bdf, 'beta')
    bdf = 1;
    pc = -1 / (2 * opts.bdf.tau * opts.bdf.beta);
else
    bdf = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute projector matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.LDL_T
    [Z, ~] = ...
        mess_column_compression(xopts.projection.Z, 'N', ...
                                xopts.projection.D, eps, 0);

    switch type
        case 'LE'

            U = 1;
            S = eqn.T;

        case 'CARE'

            if eqn.type == 'T'
                [U, S] = eig(eqn.Q);
                R = eqn.R;
            else
                [U, S] = eig(eqn.R);
                R = eqn.Q;
            end

    end

else

    [nZ, mZ] = size(xopts.projection.Z);

    if eqn.type == 'T'
        R = eye(size(eqn.B, 2));
    else
        R = eye(size(eqn.C, 1));
    end

    if xopts.projection.ortho || mZ > nZ

        Z = orth(xopts.projection.Z);
        U = 1.0;
        S = 1.0;

    else

        Z = xopts.projection.Z;
        [U, S, ~] = svd(full(Z' * Z));
        s = diag(S);
        sk = find(s > eps * s(1), 1, 'last');
        Z = Z * U(:, 1:sk) * diag(1 ./ sqrt(s(1:sk)));

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Coefficient matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type
    case 'LE'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The Lyapunov equation case
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lyapunov = true;
        C = [];
        R = [];
        B = Z' * eqn.W;
        if bdf || rosenbrock
            A = oper.mul_ApE(eqn, opts, eqn.type, ...
                             pc, eqn.type, Z, 'N');
            if bdf
                A = (opts.bdf.tau * opts.bdf.beta) * A;
                if eqn.haveUV
                    if eqn.type == 'T'
                        A = A + eqn.V * (eqn.U' * Z);
                    else
                        A = A + eqn.U * (eqn.V' * Z);
                    end
                end
            else % rosenbrock
                if opts.rosenbrock.stage == 2
                    A = (opts.rosenbrock.tau * opts.rosenbrock.gamma) * A;
                end
                if eqn.haveUV
                    if eqn.type == 'T'
                        A = A + eqn.V * (eqn.U' * Z);
                    else
                        A = A + eqn.U * (eqn.V' * Z);
                    end
                end
            end
        else
            A = oper.mul_A(eqn, opts, eqn.type, Z, 'N');
            if eqn.haveUV
                if eqn.type == 'T'
                    A = A + eqn.V * (eqn.U' * Z);
                else
                    A = A + eqn.U * (eqn.V' * Z);
                end
            end
        end
        A = Z' * A;
        if eqn.haveE
            E = Z' * oper.mul_E(eqn, opts, eqn.type, Z, 'N');
        else
            if exist('OCTAVE_VERSION', 'builtin')
                E = eye(size(A));
            else
                E = [];
            end
        end

    case 'CARE'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The Riccati equation case
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lyapunov = false;
        if eqn.type == 'T'
            opAE = 'N';
            if bdf
                B = Z' * eqn.B * sqrt(opts.bdf.tau * opts.bdf.beta);
            else
                B = Z' * eqn.B;
            end
            C = U' * eqn.C * Z;
        else
            opAE = 'T';
            if bdf
                B = Z' * eqn.C' * sqrt(opts.bdf.tau * opts.bdf.beta);
            else
                B = Z' * eqn.C';
            end
            C = U' * eqn.B' * Z;
        end
        if bdf
            A = (opts.bdf.tau * opts.bdf.beta) * ...
                (Z' * oper.mul_ApE(eqn, opts, opAE, pc, opAE, Z, 'N'));
        else
            A = oper.mul_A(eqn, opts, opAE, Z, 'N');
            if eqn.haveUV && eqn.sizeUV1
                one = 1:eqn.sizeUV1;
                if eqn.type == 'T'
                    A = A + eqn.U(:, one) * (eqn.V(:, one)' * Z);
                else
                    A = A + eqn.V(:, one) * (eqn.U(:, one)' * Z);
                end
            end
            A = Z' * A;
        end
        if eqn.haveE
            E = Z' * oper.mul_E(eqn, opts, opAE, Z, 'N');
        else
            if exist('OCTAVE_VERSION', 'builtin')
                E = eye(size(A));
            else
                E = [];
            end
        end
end

end
