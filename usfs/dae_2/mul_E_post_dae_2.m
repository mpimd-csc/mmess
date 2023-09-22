function [eqn, opts, oper] = mul_E_post_dae_2(eqn, opts, oper)
%% function post finalizes data and/or functions
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

if (not(isfield(eqn, 'Mcount'))) || not(isnumeric(eqn.Mcount))
    mess_err(opts, 'error_arguments', ['field eqn.Mcount is not defined. Did ' ...
                                       'you forget to run mul_E_pre?']);
end
if eqn.Mcount > 1
    eqn.Mcount = eqn.Mcount - 1;
else
    eqn = rmfield(eqn, 'M_');
    eqn = rmfield(eqn, 'Mcount');
end

[eqn, opts, oper] = mul_Pi_post(eqn, opts, oper);

end
