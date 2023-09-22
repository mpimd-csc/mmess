% The second order system
%
%       M x"(t) + E*x'(t) + K*x(t) = B2*u(t)
%                             y(t) = Cp*x(t) + Cv*x'(t)
%
% is implicitly transformed to the first order system
%
%                        E_f*z'(t) = A_f*z(t) + B_f*u(t)
%                             y(t) = C_f*z(t)
%
% where
%
%         | E  M|
%   E_f = | M  0|
%
%          |-K  0|
%   A_ f = | 0  M|
%
%         |  B2 |
%   B_f = |  0  |
%
%   C_f = |Cp  Cv|
%
%         | x(t)  |
%   z(t)= | x'(t) | .
%
% Matrices M, E, K are assumed to be square, symmetric and positive
% definite.
%
% The fieldnames have to end with _  to indicate that the data
% are inputdata for the algorithm.
% eqn.M_ = M
% eqn.K_ = K
% eqn.E_ = E
% eqn.B  = B_f
% eqn.C  = C_f
%
% [Benner, Kuerschner, Saak: An improved numerical method for balanced
% truncation for symmetric second-order systems, 2013]

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
