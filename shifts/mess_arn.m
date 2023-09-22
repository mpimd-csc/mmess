function [H, V, eqn, opts, oper] = mess_arn(eqn, opts, oper, opA)
%
%  Arnoldi method w.r.t. opA(A)
%
%  Calling sequence:
%
%    [H, V, eqn, opts, oper] = mess_arn(eqn, opts, oper, opA)
%
%  Input:
%
%    eqn       eqn contains data for A (Input matrix) for which
%              the Arnoldi Algorithm is supposed to run.
%    opts      struct contains parameters for the algorithm
%    oper      contains function handles with operations for A and E
%    opA       character specifies the form of opA(A)
%              opA = 'N' performs Arnoldi method
%              opA = 'I' performs inverse Arnoldi method
%
%  Output:
%
%    H         matrix H ((k+1)-x-k matrix, upper Hessenberg);
%    V         matrix V (n-x-(k+1) matrix, orthogonal columns).
%
%  Method:
%
%    The Arnoldi method produces matrices V and H such that
%
%      V(:,1) in span{b0},
%      V'*V = eye(k+1),
%      F*V(:,1:k) = V*H.
%
%    b0 = opts.shifts.b0
%
%  Remark:
%
%    This implementation does not check for (near-)breakdown!
%
%   uses operator functions size, sol_A, mul_A, sol_E, mul_E

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% Input data not completely checked!e
%% check input data
if not(isfield(opts, 'shifts')) || not(isstruct(opts.shifts))
    mess_warn(opts, 'control_data', ...
              ['shift parameter control structure missing.', ...
               'Switching to default num_desired = 25, ', ...
               'num_Ritz = 50, num_hRitz = 25.']);
    opts.shifts.num_desired = 25;
    opts.shifts.num_Ritz = 50;
    opts.shifts.num_hRitz = 25;
else
    if not(isfield(opts.shifts, 'num_desired')) || ...
            not(isnumeric(opts.shifts.num_desired))
        mess_warn(opts, 'control_data', ...
                  ['Missing or corrupted opts.shifts.num_desired ', ...
                   'field. Switching to default: 25']);
        opts.shifts.num_desired = 25;
    end
    if strcmp(opts.shifts.method, 'heur') && ...
       (not(isfield(opts.shifts, 'num_Ritz')) || ...
        not(isnumeric(opts.shifts.num_Ritz)))
        mess_warn(opts, 'control_data', ...
                  ['Missing or Corrupted opts.shifts.num_Ritz ', ...
                   'field. Switching to default: 50']);
        opts.shifts.num_Ritz = 50;
    end
    if strcmp(opts.shifts.method, 'heur') && ...
       (not(isfield(opts.shifts, 'num_hRitz')) || ...
        not(isnumeric(opts.shifts.num_hRitz)))
        mess_warn(opts, 'control_data', ...
                  ['Missing or Corrupted opts.shifts.num_hRitz ', ...
                   'field. Switching to default: 25']);
        opts.shifts.num_hRitz = 25;
    end
end

if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');
if not(result)
    mess_err(opts, 'control_data', ...
             'system data is not completely defined or corrupted');
end
if isfield(eqn, 'haveUV') && eqn.haveUV
    if not(isfield(eqn, 'U')) || isempty(eqn.U) || ...
            not(isfield(eqn, 'V')) || isempty(eqn.V) || ...
            not((size(eqn.U, 1)) == size(eqn.V, 1) && ...
                size(eqn.U, 2) == size(eqn.V, 2))
        mess_err(opts, 'SMW', ...
                 'Inappropriate data of low-rank updated ', ...
                 'operator (eqn.U and eqn.V)');
    end
end
if not(isfield(opts, 'rosenbrock'))
    opts.rosenbrock = [];
end
if isstruct(opts.rosenbrock) && isfield(opts.rosenbrock, 'tau')
    rosenbrock = true;
    pc = -1 / (2 * opts.rosenbrock.tau);
else
    rosenbrock = false;
end
if not(isfield(opts, 'bdf'))
    opts.bdf = [];
end
if isstruct(opts.bdf) && isfield(opts.bdf, 'tau') && ...
        isfield(opts.bdf, 'beta')
    bdf = true;
    pc = -1 / (2 * opts.bdf.tau * opts.bdf.beta);
else
    bdf = false;
end

%% check input Parameters
if not(ischar(opA))
    mess_err(opts, 'error_arguments', 'opA is not a char');
end

opA = upper(opA);
if not(opA == 'N' || opA == 'I')
    mess_err(opts, 'error_arguments', 'opA is not ''N'' or ''I''');
end

% returns order of A or states of A, A is supposed to be square
n = oper.size(eqn, opts);

if opA == 'I'
    k = opts.shifts.num_hRitz;
else
    k = opts.shifts.num_Ritz;
end

if k >= n - 1
    mess_err(opts, 'error_arguments', ...
             'k must be smaller than the order of A!');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_E_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_A_pre(eqn, opts, oper);

if bdf || rosenbrock
    [eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);
end
%% initialize data
if not(isfield(opts.shifts, 'b0')) || isempty(opts.shifts.b0)
    b0 = ones(n, 1);
else
    b0 = opts.shifts.b0;
end

V = zeros(length(b0), k + 1);
H = zeros(k + 1, k);

V(:, 1) = (1.0 / norm(b0)) * b0;

beta = 0;
%% perform Arnoldi method
for j = 1:k

    if j > 1
        V(:, j) = (1.0 / beta) * w;
    end

    % no eqn.type cases needed, eigenvalues are the same for transposed
    % operator
    if opA == 'I' % Perform inverse Arnoldi
        if isfield(eqn, 'haveE') && eqn.haveE

            EV = oper.mul_E(eqn, opts, 'N', V(:, j), 'N');

            if isfield(eqn, 'haveUV') && eqn.haveUV
                RHS = [EV eqn.U];
                if bdf
                    AB = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                                      -2 * pc * RHS, 'N');
                elseif rosenbrock
                    AB = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                                      RHS, 'N');
                else
                    AB = oper.sol_A(eqn, opts, 'N', RHS, 'N');
                end
                AV = AB(:, 1);
                AU = AB(:, 2:end);
                Im = speye(size(eqn.U, 2));
                w = AV - AU * ((Im + eqn.V' * AU) \ (eqn.V' * AV));
            else % no rosenbrock case here as that always has UV
                if bdf
                    w = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                                     -2 * pc * EV, 'N');
                else
                    w = oper.sol_A(eqn, opts, 'N', EV, 'N');
                end
            end
        else
            if isfield(eqn, 'haveUV') && eqn.haveUV
                RHS  = [V(:, j) eqn.U];
                if bdf
                    AB = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                                      -2 * pc * RHS, 'N');
                elseif rosenbrock
                    AB = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                                      RHS, 'N');
                else
                    AB = oper.sol_A(eqn, opts, 'N', RHS, 'N');
                end
                AV = AB(:, 1);
                AU = AB(:, 2:end);
                Im = speye(size(eqn.U, 2));
                w = AV - AU * ((Im + eqn.V' * AU) \ (eqn.V' * AV));
            else % no rosenbrock case here as that always has UV
                if bdf
                    w = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                                     -2 * pc * V(:, j), 'N');
                else
                    w = oper.sol_A(eqn, opts, 'N', V(:, j), 'N');
                end
            end
        end
    else % opA = 'N' Perform standard Arnoldi
        if isfield(eqn, 'haveE') && eqn.haveE
            if isfield(eqn, 'haveUV') && eqn.haveUV
                UVtV =  eqn.U * (eqn.V' * V(:, j));
                if bdf || rosenbrock
                    AV = oper.mul_ApE(eqn, opts, 'N', pc, 'N', ...
                                      V(:, j), 'N');
                else
                    AV = oper.mul_A(eqn, opts, 'N', V(:, j), 'N');
                end

                if bdf
                    w = oper.sol_E(eqn, opts, 'N', ...
                                   (-2 / pc) * AV + UVtV, 'N');
                else % rosenbrock and default are now doing the same
                    w = oper.sol_E(eqn, opts, 'N', ...
                                   AV + UVtV, 'N');
                end
            else % no rosenbrock case here as that always has UV
                if bdf
                    ApEV = oper.mul_ApE(eqn, opts, 'N', pc, 'N', ...
                                        V(:, j), 'N');
                    w = oper.sol_E(eqn, opts, 'N', (-2 / pc) * ...
                                   ApEV, 'N');
                else
                    w = oper.sol_E(eqn, opts, 'N', ...
                                   oper.mul_A(eqn, opts, 'N', ...
                                              V(:, j), 'N'), 'N');
                end
            end
        else
            if isfield(eqn, 'haveUV') && eqn.haveUV
                UVtV =  eqn.U * (eqn.V' * V(:, j));
                if bdf
                    w = (-2 / pc) * ...
                        oper.mul_ApE(eqn, opts, 'N', pc, 'N', ...
                                     V(:, j), 'N');
                elseif rosenbrock
                    w = oper.mul_ApE(eqn, opts, 'N', pc, 'N', ...
                                     V(:, j), 'N');
                else
                    w = oper.mul_A(eqn, opts, 'N', V(:, j), 'N');
                end
                w = w + UVtV;
            else % no rosenbrock case here as that always has UV
                if bdf
                    w = (opts.bdf.tau * opts.bdf.beta) * ...
                        oper.mul_ApE(eqn, opts, 'N', pc, 'N',  ...
                                     V(:, j), 'N');
                else
                    w = oper.mul_A(eqn, opts, 'N', V(:, j), 'N');
                end
            end
        end
    end

    for k = 1:2 % repeated MGS
        for i = 1:j
            g = V(:, i)' * w;
            H(i, j) = H(i, j) + g;
            w = w -   V(:, i) * g;
        end
    end
    beta = norm(w);
    H(j + 1, j) = beta;

end

V(:, k + 1) = (1.0 / beta) * w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_E_post(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_A_post(eqn, opts, oper);

if bdf || rosenbrock
    [eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);
end

end
