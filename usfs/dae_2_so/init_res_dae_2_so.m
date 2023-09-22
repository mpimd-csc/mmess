function [W, res0, eqn, opts, oper] = init_res_dae_2_so(eqn, opts, oper, W, T)
%% function init_res initializes the low-rank residual W and res0
% function [ W, res0, eqn, opts, oper ] = ...
%      init_res_dae_2_so( eqn, opts, oper, W, T)
%
%   Input/Output:
%
%   eqn        structure containing data for G or B or C
%   opts       structure containing parameters for the algorithm
%   oper       struct contains function handles for operation with A and E
%   W          right hand side matrix
%   T          matrix such that the residual is W*T*W'
%              (optional, defaults to identity
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
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check data
for mat = 'MEKG'
    if not(isfield(eqn, sprintf('%c_', mat))) || not(eval(sprintf('isnumeric(eqn.%c_)', mat)))
        mess_err(opts, 'error_arguments', 'field eqn.%c_ is not defined', mat);
    end
end

nv = size(eqn.M_, 1);
np = size(eqn.G_, 1);

if not(isfield(eqn, 'type'))
    eqn.type = 'N';
    mess_warn(opts, 'equation_type', ['Unable to determine type of equation.'...
                                      'Falling back to type ''N''']);
end

if (not(isnumeric(W))) || (not(ismatrix(W)))
    mess_err(opts, 'error_arguments', 'W has to be a matrix');
end

%% compute low-rank residual
Wtemp = zeros(2 * nv + np, size(W, 2));
Wtemp(1:size(W, 1), :) = W;

% if eqn.type=='N'
%    S = [speye(nv,nv),sparse(nv,nv),sparse(nv,np);
%          sparse(nv,nv),eqn.M_,eqn.G_';
%          sparse(np,nv),eqn.G_,sparse(np,np)];
%     X = full( S \ Wtemp);
%     W = [X(1:nv,:);eqn.M_*X(nv+1:2*nv,:)];
% else
%     S = [speye(nv,nv),sparse(nv,nv),sparse(nv,np);
%          sparse(nv,nv),eqn.M_',eqn.G_';
%          sparse(np,nv),eqn.G_,sparse(np,np)];
%     %X = full( S \ [W; sparse(np,size(W, 2))]);
%     X = full( S \ Wtemp);
%     W = [X(1:nv,:);eqn.M_'*X(nv+1:2*nv,:)];
% end

S = [speye(nv, nv), sparse(nv, nv), sparse(nv, np)
     sparse(nv, nv), eqn.M_, eqn.G_'
     sparse(np, nv), eqn.G_, sparse(np, np)];
X = full(S \ Wtemp);
W = [X(1:nv, :); eqn.M_ * X(nv + 1:2 * nv, :)];

%% compute res0
if not(exist('T', 'var')) && opts.LDL_T
    % this means we only use init_res for potential projection
    return
end

if isfield(opts, 'nm') && isfield(opts.nm, 'res0')
    res0 = opts.nm.res0;
else
    if opts.LDL_T
        if opts.norm == 2
            res0 = max(abs(eig(W' * W * T)));
        else
            res0 = norm(eig(W' * W * T), 'fro');
        end
    else
        res0 = norm(W' * W, opts.norm);
    end
end

end
