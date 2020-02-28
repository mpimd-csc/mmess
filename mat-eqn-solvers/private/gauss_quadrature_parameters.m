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

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009-2020
%

% Golub-Welsch algorithm on the interval [-1, 1]
N = (order - 1) / 2;
j = 1:N;
beta = j ./ sqrt(4*j.^2 - 1);

J = diag(beta(1:end-1), -1) + diag(beta(1:end-1), 1);
[V, nodes] = eig(J, 'vector');

[nodes, ind] = sort(nodes, 'ascend');
weights = V(1, ind).^2 * 2;

% Conversion from [-1, 1] to [0, h]:
weights = weights' * h/2;
nodes = h/2*(nodes + 1);

end

