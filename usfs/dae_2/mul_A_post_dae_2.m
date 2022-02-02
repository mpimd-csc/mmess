function [eqn,opts,oper] = mul_A_post_dae_2(eqn,opts,oper)
%  MUL_A_POST_DAE_2 clears the hidden manifold projector used in
%  mul_A_dae_2.
%
% Input/Output:
%  eqn, opts, oper  the usual equation, options and function handle
%                   structures; passed trough to enable clearing of
%                   helper data inside them.
%

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

[eqn,opts,oper] = mul_Pi_post(eqn,opts,oper);