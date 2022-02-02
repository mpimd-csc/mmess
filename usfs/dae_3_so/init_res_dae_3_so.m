function [ W, res0, eqn, opts, oper ] = init_res_dae_3_so( eqn, opts, oper, RHS)
%% function init_res initializes the low rank residual W and res0
% function [ W, res0, eqn, opts, oper ] = init_res_dae_3_so( eqn, opts, oper, RHS)
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
%   uses no other dae_3_so function

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check data
for mat='MEKG'
    if(not(isfield(eqn, sprintf('%c_',mat))) || not(eval(sprintf('isnumeric(eqn.%c_)',mat))))
        error('MESS:error_arguments', 'field eqn.%c_ is not defined',mat);
    end
end

nv = size(eqn.M_,1);
np = size(eqn.G_,1);

if not(isfield(eqn,'type'))
  eqn.type='N';
  warning('MESS:equation_type',['Unable to determine type of equation.'...
    'Falling back to type ''N''']);
end

if (not(isnumeric(RHS))) || (not(ismatrix(RHS)))
    error('MESS:error_arguments','RHS has to be a matrix');
end

%% compute low rank residual
RHStemp1 = zeros(nv+np,size(RHS,2));
RHStemp2 = zeros(nv+np,size(RHS,2));
RHStemp1(1:nv,:) = RHS(1:nv,:);
RHStemp2(1:nv,:) = RHS(nv+1:2*nv,:);

S = [eqn.M_,eqn.G_';
     eqn.G_,sparse(np,np)];
X1 = full( S \ RHStemp1);
X2 = full( S \ RHStemp2);
W = [eqn.M_*X1(1:nv,:);eqn.M_*X2(1:nv,:)];



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
        res0 = norm(RHS' * RHS, opts.norm);
    end
end

end

