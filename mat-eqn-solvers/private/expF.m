function [L, D, eqn, opts, oper] = expF(eqn, opts, oper, h, IQL, IQD, L0, D0, t0)
% Solve the affine problem arising from the split DRE.
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
if nargin < 9
    t0 = 0;
end

[out, eqn, opts, oper] = mess_exp_action(eqn, opts, oper, h, L0, t0);

L = [out.Z, IQL];
D = blkdiag(D0, IQD);

[L, D] = mess_column_compression(L, 'N', D, opts.splitting.trunc_tol, ...
    opts.splitting.trunc_info);
end