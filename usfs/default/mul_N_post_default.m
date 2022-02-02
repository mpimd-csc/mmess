function [eqn, opts, oper] = mul_N_post_default(eqn, opts, oper)
% Transforms eqn.N_ back into a matrix (if it was given as such) and
% if it the last call of mul_N_post (eqn.Ncount = 1)
%
% Input/Output:
%    eqn    struct contains data for equations
%
%    opts   struct contains parameters for the algorithm
%
%    oper   struct contains function handles for operation with N
%
%
% input        eqn.N_           (as cell or matrix)
%              eqn.originalN   (saves matrix version for post_N)
%              eqn.Ncount      (Function calls of mul_N_pre)
%
% output       eqn.N_           (as cell or matrix)
%              eqn.Ncount      (should be 1)

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


% checks if mul_E_pre was initialized
if(not(isfield(eqn, 'Ncount'))) || not(isnumeric(eqn.Ncount))
    error('MESS:error_arguments', ['field eqn.Ncount is not defined. Did ' ...
        'you forget to run mul_E_pre?']);
end

% checks Ncount and decides output as cell or matrix
if eqn.Ncount > 1
    eqn.Ncount = eqn.Ncount - 1;
else
    if not(isfield(eqn, 'originalN'))||isempty(eqn.originalN)
        eqn.N_ = eqn.N_;
    else
        eqn.N_ = eqn.originalN;
    end
end

end