% Function Handles for structured index-1 differential-algebraic equations
% of second order, see e.g., [1].
%
% Differential-Algebraic System
% | M11 0 |          | E11 0 |         | K11 K12 |        | B1 |
% |       | p''(t) + |       | p'(t) + |         | p(t) + |    | u(t),
% |  0  0 |          |  0  0 |         | K21 K22 |        | B2 |
%                                                                  (1)
%            y(t) = | C1 C2 | p(t)
%
% Attention, the matrices M11, E11 and A22 need to be invertible.
%
% Matrices in the structure are
% eqn.M_
% eqn.E_
% eqn.K_
% eqn.B
% eqn.C
%
% in exactly the form above.
%
% The size of the square matrices M11 and E11 coincides and is stored in
% eqn.nd.
%
% Implicitly the system is lifted to first order form
%
% | E11 M11 |         | -K   0  |        | B |
% |         | x'(t) = |         | x(t) + |   | u(t).
% | M11  0  |         | 0   M11 |        | 0 |
%                                                                (2)
%              y(t) = | C 0 | x(t) + D u(t),
%
% where K = K11 - K12 * K22 \ K21,   B = B1 - K12 * K22 \B2,
% C = C1 - C2 * K22 \ K21 and D = C2 * K22 \ B2.
%
% Note that eqn.B and eqn.C are overwritten by their corresponding
% representations on the 2*eqn.nd dimensional hidden manifold, i.e. in the
% first order ODE realization of the system.
%
% References
%
%   [1] P. Benner, J. Saak, M. M. Uddin, Structure preserving model order
%       reduction of large sparse second-order index-1 systems and
%       application to a mechatronics model,
%       Math. Comput. Model. Dyn. Syst. 22 (6) (2016) 509–523.
%       https://doi.org/10.1080/13873954.2016.1218347.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
