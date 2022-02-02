% Function Handles for structured index-1 differential-algebraic equations,
% e.g., power systems examples from
% https://morwiki.mpi-magdeburg.mpg.de/morwiki/index.php/Power_system_examples
% 
% Differential-Algebraic System
% | E11 0 |         | A11 A12 |        | B1 |
% |       | x'(t) = |         | x(t) + |    | u(t),
% |  0  0 |         | A21 A22 |        | B2 |
% 
%            y(t) = | C1 C2 | x(t)
% 
% Attention, the matrices E11 and A22 need to be invertible.
% The fieldnames have to end with _ to indicate that the Data are inputdata
% for the Algorithm:
% eqn.A_
% eqn.E_
% eqn.B
% eqn.C
% 
% Note that eqn.B and eqn.C are overwritten by their corresponding 
% representations on the hidden manifold, i.e. in the ODE realization of 
% the system.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
