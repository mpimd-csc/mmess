function [ eqn, opts, oper ] = init_res_pre_dae_2( eqn, opts, oper )
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
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


[eqn, opts, oper] = mul_Pi_pre(eqn,opts,oper);
end

