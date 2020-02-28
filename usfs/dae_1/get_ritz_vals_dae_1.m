function [rw, Hp, Hm, Vp, Vm] = get_ritz_vals_dae_1(eqn, opts, oper, U, W, p_old)
%  This function ensures that W is not empty if the shift
%  method is projection. Otherwise, it checks opts.shifts.b0. It 
%  should be a vector of the same size as eqn.A but oper.size gives eqn.st.
%
%  MMESS (Björn Baran, March 2018)

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

% Input data not completely checked!
if(not(isfield(eqn,'A_'))) || not(isnumeric(eqn.A_))
    error('MESS:error_arguments','field eqn.A_ is not defined');
end
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');
if not(result)
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end
n = oper.size(eqn, opts);

%%

if isfield(opts.shifts, 'method') && ...
        strcmp(opts.shifts.method, 'projection')
    if isempty(W)
        % first shifts are computed with U = eqn.G and W = A * eqn.G
        W = oper.mul_A(eqn, opts, eqn.type, U, 'N');
    end
    rw = mess_projection_shifts(eqn, opts, oper, U, ...
        W, p_old);
else
    if (not(isfield(opts.shifts, 'b0')) || isempty(opts.shifts.b0))
        opts.shifts.b0 = ones(n,1);
    else
        if length(opts.shifts.b0) ~= n
            warning('MESS:b0',...
                'b0 has the wrong length. Switching to default.');
            opts.shifts.b0 = ones(n,1);
        end
    end
    [rw, Hp, Hm, Vp, Vm] = mess_get_ritz_vals(eqn, opts, oper);
end
