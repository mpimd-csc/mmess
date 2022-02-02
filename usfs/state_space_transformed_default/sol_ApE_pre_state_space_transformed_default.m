function [eqn, opts, oper] = ...
    sol_ApE_pre_state_space_transformed_default(eqn, opts, oper)
%% function [eqn, opts, oper] = ...
%     sol_ApE_pre_state_space_transformed_default(eqn, opts, oper)
%
% function pre initializes data and/or functions
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
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


if isfield(eqn, 'EL') && isfield(eqn, 'EU')
    if isfield(eqn, 'Ecount')
        eqn.Ecount = eqn.Ecount + 1;
    else
        eqn.Ecount = 2;
    end
else
    [eqn.EL, eqn.EU] = lu(eqn.E_);
    eqn.Ecount       = 1;
end
