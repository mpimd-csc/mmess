function [ W, res0, eqn, opts, oper ] = init_res_dae_2( eqn, opts, oper, RHS)
%% function init_res initializes the low rank residual W and res0
% function [ W, res0, eqn, opts, oper ] = init_res_dae_2( eqn, opts, oper, RHS)
%
%   Input/Output:
%
%   eqn        structure containing data for G or B or C
%   opts       structure containing parameters for the algorithm
%   oper       struct contains function handles for operation with A and E
%   RHS        right hand side matrix
%
%   Outputs:
%
%   W          matrix given by ADI to compute residuum
%   res0       initial residuum norm
%
%   uses no other dae_2 function

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check data
if not(isfield(eqn,'A_')) || not(isnumeric(eqn.A_))
    error('MESS:equation_data',...
      'Empty or Corrupted field A detected in equation structure.');
end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end
if not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_))
    error('MESS:equation_data',...
      'Empty or Corrupted field E detected in equation structure.');
end
if (not(isnumeric(RHS))) || (not(ismatrix(RHS)))
    error('MESS:error_arguments','RHS has to ba a matrix');
end
if (eqn.st ~= size(RHS, 1))
    error('MESS:error_arguments','eqn.st differs with number of rows of RHS');
end

%% compute low rank residual

W = mul_Pi(eqn,'N',RHS,'N');

%% compute res0
if isfield(opts, 'nm') && isfield(opts.nm, 'res0')
    res0 = opts.nm.res0;
else
    if opts.LDL_T
        if opts.norm == 2
            res0 = max(abs(eig(RHS' * RHS * diag(eqn.S_diag))));
        else
            res0 = norm(eig(RHS' * RHS * diag(eqn.S_diag)), 'fro');
        end
    else
        res0 = norm(full(RHS' * RHS), opts.norm); %opts.norm == 2 needs dense matrix
    end
end

end

