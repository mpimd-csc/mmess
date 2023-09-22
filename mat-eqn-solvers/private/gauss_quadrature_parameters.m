function [nodes, weights] = gauss_quadrature_parameters(h, order)
% Return nodes (x) and weights (w) of the Gauss-Legendre quadrature rule
% of order 2N+1, i.e. the parameter order must be odd >= 3.
% The interval is [0, h] rather than [-1, 1]. Uses the Golub-Welsch [1]
% algorithm.
%
% [1] Gene H. Golub and John H. Welsch, Calculation of Gauss quadrature
% rules, Math. Comp. 23 (1969), pp. 221-230. DOI:
% https://doi.org/10.1090/S0025-5718-69-99647-1

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

% Golub-Welsch algorithm on the interval [-1, 1]
N = (order - 1) / 2;
k = 1:N;
beta = k ./ sqrt(4 * k.^2 - 1);

J = diag(beta(1:end - 1), -1) + diag(beta(1:end - 1), 1);
[V, nodes] = eig(J, 'vector');

[nodes, ind] = sort(nodes, 'ascend');
weights = V(1, ind).^2 * 2;

% Conversion from [-1, 1] to [0, h]:
weights = weights' * h / 2;
nodes = h / 2 * (nodes + 1);

end
