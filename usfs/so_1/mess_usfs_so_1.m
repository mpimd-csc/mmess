function mess_usfs_so_1
%% function mess_usfs_so_1
% Second Order System
%
%        M*x"(t) + E*x'(t) + K*x(t) = B2*u(t)
%                              y(t) = Cp*x(t) + Cv*x'(t)
%
%  is transformed to the first order system
%
%                          E_f*x_f' = A_f*x_f + B_f*u
%                              y(t) = C_f*x_f
%
%  where
%
%           |-K  0 |
%     E_f = | 0  M | ,
%
%           | 0 -K |
%     A_f = |-K -E |,
%
%           | 0  |
%     B_f = | B2 |,
%
%     C_f = |Cp Cv|
%
%           | x |
%     x_f = | x'|.
%
% So we have:
%
% |-K 0||x'|= |0  -K||x | + |0|
% |0  M||x"|  |-K -E||x'|   |B|u
%
%
% Attention the Matrix M E K are symmetric and quadratic.
% K is a fullrank Matrix.
% The fieldnames have to end with _  to indicate that the Data
% are inputdata for the Algorithm.
% eqn.M_ = M
% eqn.K_ = K
% eqn.E_ = E
% eqn.B  = B_f
% eqn.C  = C_f

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
