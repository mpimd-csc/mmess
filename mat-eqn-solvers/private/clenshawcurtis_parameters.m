function [x, w] = clenshawcurtis_parameters(a, b, N)
% Return nodes (x) and weights (w) of a Clenshaw-Curtis quadrature rule
% with N+1 points on the interval [a, b]. N must be even.
% The rule is of order N (error prop. to (a-b)^{-N}), but typically 
% performs better. The nodes are interleaved in the sense that if 
% N and N2 = 2N have the corresponding nodes x and x2, then 
% x2(1:2:end) = x; This can be used for error estimation.
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

h = b-a; % interval length

% Column 1 used for the weights, column 2 for the nodes, so we don't have
% to call cos() at all (by setting frequency 1 to 1 and the rest to 0)
c = zeros(N+1, 2);
c(1:2:N+1,1) = (2./[1 1-(2:2:N).^2 ])';
c(2,2) = 1;

xi = real(ifft([c(1:N+1,:); c(N:-1:2,:)])); % symmetrize

w = h*xi(1:N+1,1);
w(1) = w(1) / 2;
w(end) = w(end) / 2;

x = ((b+a)/2 + N*h/2*xi(1:N+1,2));

end
