function s = outerfrobnormdiff_LDLT(L1, D1, L2, D2)
% Compute the Frobenius norm of L1 D1 L1' - L2 D2 L2' 
% for symmetric D1, D2
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

    [L, D] = mess_column_compression([L1, L2], 'N', blkdiag(D1, -D2), ...
        eps, 0);
    % Exploit the invariance under cyclic permutations property of the trace
    % to compute norm(L*D*L','fro')
    s = sqrt(trace((L'*L*D)^2)); 

    % Expensive version
%     s = norm(L1*D1*L1' - L2*D2*L2', 'fro');

end
