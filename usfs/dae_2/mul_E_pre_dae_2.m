function [eqn, opts, oper] = mul_E_pre_dae_2(eqn, opts, oper)
%% function pre initializes data and/or functions
%
% Input:
%    eqn    struct contains data for equations
%
%    opts   struct contains parameters for the algorithm
%
%    oper   struct contains function handles for operation with A
%
% Output:
% eqn
% opts
% oper

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

alpha = -1 / 50;
if isfield(eqn, 'manifold_dim') && isnumeric(eqn.manifold_dim)
    one = 1:eqn.manifold_dim;
else
    mess_err(opts, 'wrong_arguments', ...
             'missing or corrupted field st detected');
end
if not(isfield(eqn, 'M_'))
    if not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_)) || ...
            not(isfield(eqn, 'A_')) || not(isnumeric(eqn.A_))
        mess_err(opts, 'error_arguments', ...
                 'field eqn.E_ or eqn.A_ is not defined or corrupted');
    end
    eqn.M_ =           alpha * eqn.A_;
    eqn.M_(one, one) = eqn.E_(one, one);
    eqn.Mcount =       1;
else
    if not(isfield(eqn, 'Mcount')) || not(isnumeric(eqn.Mcount))
        mess_err(opts, 'error_arguments', ...
                 'field eqn.Mcount is not defined. Did ', ...
                 'you forget to run mul_E_pre?');
    end
    eqn.Mcount = eqn.Mcount + 1;
end

[eqn, opts, oper] = mul_Pi_pre(eqn, opts, oper);
end
