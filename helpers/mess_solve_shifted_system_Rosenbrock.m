function [ V, eqn, opts, oper ]=mess_solve_shifted_system_Rosenbrock(eqn, opts, oper, pc, W)
% Solves (Ã + p*E)V = W for V, Ã = A - 1/(2*tau)*E - UV^T
%  (Rosenbrock Scheme)
%
%  Solves (Ã + p*E)V = W for V, Ã = A - 1/(2*tau)*E - UV^T if eqn.type == 'N'
%  Solves (Ã + p*E)^T*V = W for V, Ã = A - 1/(2*tau)*E - UV^T if eqn.type == 'T'
%   (Rosenbrock Scheme)
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
% Output:
%  V         solution of the shifted system
%
%  eqn       structure containing equation data
%
%  opts      structure containing parameters for the algorithm
%
%  oper      contains function handles with operations for A and E

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Check input

%% Initialize data
k = size(W, 2);
if eqn.haveUV
    m = size(eqn.U, 2);
end

%% Manipulate shift for Rosenbrock scheme
if opts.rosenbrock.stage == 1
    pc = pc - 1 / (opts.rosenbrock.tau * 2);
else % p = 2
    taugamma = (opts.rosenbrock.tau * opts.rosenbrock.gamma);
    pc = (pc - 0.5) / taugamma;
end

%% preprocess shifted solver
[eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);

%% solve shifted system
if opts.rosenbrock.stage == 1  % 1st order Rosenbrock
    if eqn.haveUV %Perform Sherman-Morrison-Woodbury-trick
        if eqn.type == 'T'
            V = oper.sol_ApE(eqn, opts,eqn.type,pc,eqn.type,[W eqn.V],'N');
            SMW = V(:,k+1:end);
            V=V(:,1:k);
            V=V-SMW*((eye(m)+eqn.U'*SMW)\(eqn.U'*V));
        else
            V = oper.sol_ApE(eqn, opts,eqn.type,pc,eqn.type,[W eqn.U],'N');
            SMW = V(:,k+1:end);
            V=V(:,1:k);
            V=V-SMW*((eye(m)+eqn.V'*SMW)\(eqn.V'*V));

        end
    else
        V = oper.sol_ApE(eqn, opts,eqn.type,pc,eqn.type,W,'N');
    end
else % p = 2, 2nd order Rosenbrock
    if eqn.haveUV %Perform Sherman-Morrison-Woodbury-trick
        if eqn.type == 'T'
            V = oper.sol_ApE(eqn, opts,eqn.type,pc,eqn.type,[W eqn.V] / taugamma,'N');
            SMW = V(:,k+1:end);
            V=V(:,1:k);
            V=V-SMW*((eye(m)+eqn.U'*SMW)\(eqn.U'*V));
        else
            V = oper.sol_ApE(eqn, opts,eqn.type,pc,eqn.type,[W eqn.U] / taugamma,'N');
            SMW = V(:,k+1:end);
            V=V(:,1:k);
            V=V-SMW*((eye(m)+eqn.V'*SMW)\(eqn.V'*V));

        end
    else
        V = oper.sol_ApE(eqn, opts,eqn.type,pc,eqn.type, W / taugamma,'N');
    end
end

%% postprocess shifted solver
[eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);
