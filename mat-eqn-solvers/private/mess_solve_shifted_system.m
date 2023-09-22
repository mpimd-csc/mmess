function [X, eqn, opts, oper] = ...
    mess_solve_shifted_system(eqn, opts, oper, pc, W, mode)
% Solves (Ã + p*E)X = W for X, Ã = A or Ã = A + UV^T
%
%  Solves (Ã + p*E)X = W for X, Ã = A or Ã = A + UV^T if eqn.type == 'N'
%  Solves (Ã + p*E)^T*X = W for X, Ã = A or Ã = A + UV^T if eqn.type == 'T'
%
%
% Input:
%  eqn       structure containing equation data
%
%  opts      structure containing parameters for the algorithm
%
%  oper      contains function handles with operations for A and E
%
%  pc        contains shift parameter p
%
%  W         contains right hand side
%
%  mode      decides if we are running plain solvers or as part of
%            a Rosenbrock or BDF method (optional, default: empty)
%
% Output:
%  X         solution of the shifted system
%
%  eqn       structure containing equation data
%
%  opts      structure containing parameters for the algorithm
%
%  oper      contains function handles with operations for A and E

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check input
haveUV =  isfield(eqn, 'haveUV') && eqn.haveUV;
%% Initialize data
k = size(W, 2);
Wcols = 1:k;

if haveUV
    m = size(eqn.U, 2);
    UVcols = (k + 1):(k + m);
    switch eqn.type
        case 'N'
            RHS = [W, eqn.U];
            V = eqn.V;
        case 'T'
            RHS = [W, eqn.V];
            V = eqn.U;
    end
else
    RHS = W;
end

if not(exist('mode', 'var')) || isempty(mode)
    mode = 'default';
end

switch lower(mode)
    case 'bdf'
        tau_beta = opts.bdf.tau * opts.bdf.beta;
        pc = (pc - 0.5) / tau_beta;
        RHS = RHS / tau_beta;
    case 'rosenbrock'
        if opts.rosenbrock.stage == 1
            pc = pc - 1 / (opts.rosenbrock.tau * 2);
        else % p = 2
            tau_gamma = (opts.rosenbrock.tau * opts.rosenbrock.gamma);
            pc = (pc - 0.5) / tau_gamma;
            RHS = RHS / tau_gamma;
        end
end

%% preprocess shifted solver
[eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);

%% solve shifted system
X = oper.sol_ApE(eqn, opts, eqn.type, pc, eqn.type, RHS, 'N');
if haveUV
    % Perform Sherman-Morrison-Woodbury-trick
    SMW = X(:, UVcols);
    X = X(:, Wcols);
    X = X - SMW * ((eye(m) + V' * SMW) \ (V' * X));
end

%% postprocess shifted solver
[eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);
end
