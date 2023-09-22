function [x, w] = clenshawcurtis_parameters(a, b, N)
% Return nodes (x) and weights (w) of a Clenshaw-Curtis quadrature rule
% with N+1 points on the interval [a, b]. N must be even.
% The rule is of order N (error prop. to (a-b)^{-N}), but typically
% performs better. The nodes are interleaved in the sense that if
% N and N2 = 2N have the corresponding nodes x and x2, then
% x2(1:2:end) = x; This can be used for error estimation.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

h = b - a; % interval length

% Column 1 used for the weights, column 2 for the nodes, so we don't have
% to call cos() at all (by setting frequency 1 to 1 and the rest to 0)
c = zeros(N + 1, 2);
c(1:2:N + 1, 1) = (2 ./ [1 1 - (2:2:N).^2])';
c(2, 2) = 1;

xi = real(ifft([c(1:N + 1, :); c(N:-1:2, :)])); % symmetrize

w = h * xi(1:N + 1, 1);
w(1) = w(1) / 2;
w(end) = w(end) / 2;

x = ((b + a) / 2 + N * h / 2 * xi(1:N + 1, 2));

end
