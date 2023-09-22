function [eqn, opts, oper] = ...
    mul_ApE_post_state_space_transformed_default(eqn, opts, oper)
%% function [eqn, opts, oper] = ...
%     mul_ApE_post_state_space_transformed_default(eqn, opts, oper)
%
% function post finalizes data and/or functions
%
% Input
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   oper            struct contains function handles for operation
%                   with A and E
%
% Output
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   oper            struct contains function handles for operation
%                   with A and E

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

mess_assert(opts, isfield(eqn, 'Ecount'), ...
            'error_arguments', ...
            'field eqn.Ecount is not defined.');

if eqn.Ecount > 1
    eqn.Ecount = eqn.Ecount - 1;
else
    eqn = rmfield(eqn, 'Ecount');
    eqn = rmfield(eqn, 'EL');
    eqn = rmfield(eqn, 'EU');
end
