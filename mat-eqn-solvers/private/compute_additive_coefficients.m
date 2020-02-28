function [gamma] = compute_additive_coefficients(order, symmetric)
% Compute order conditions for the asymmetric and symmetric additive
% splitting schemes
% 
% The order must be even for the symmetric schemes.
% 
% Use symmetric = true for the symmetric version, symmetric = false for the
% asymmetric. 
% 

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


if symmetric
    s = order / 2;

    A = zeros(s,s);
    A(1,:) = 1;
    for k = 1:s-1
        A(k+1, :) = (1:s).^(-2*k);
    end

    b = [1/2; zeros(s-1, 1)];

    gamma = A\b;
else
    s = order;

    A = zeros(s,s);
    A(1,:) = 1;
    for k = 1:s-1
        A(k+1, :) = (1:s).^(-k);
    end

    b = [1; zeros(s-1, 1)];

    gamma = A\b;
end