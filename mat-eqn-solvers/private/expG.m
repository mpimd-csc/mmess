function [L, D] = expG(eqn, opts, oper, h, L0, D0, t0)
% Solve the nonlinear problem arising from the split DRE
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
%               2009-2019
%
L = L0;
I = eye(size(D0));

if not(eqn.LTV) % Autonomous case
    if eqn.type == 'T'
        D = (I + h*D0*L0'*(eqn.B*(eqn.Rinv*(eqn.B'*L0))))\D0;
    elseif eqn.type == 'N'
        D = (I + h*D0*L0'*(eqn.C'*(eqn.Rinv*(eqn.C*L0))))\D0;
    end
else 
    % Time-varying case. Only difference is that instead of h*B*Rinv*B' we 
    % have int_{t}^{t+h} {B(s)RinvB(s)' ds} we approximate this as 
    % LB*DB*LB' by quadrature. Note that the factor h is included in the 
    % approximation. Could probably get away with lower order...
    [xj, wj] = exact_quadrature_parameters(h, 29);
    xj = xj + t0; % Shift from [0, h] t0 [t0, t0+h]

    % Need to evaluate eqn.B or eqn.C' at all xj(j)
    LB = cell(1, length(xj));
    for j = 1:length(xj)
       [eqn, opts, oper] = ...
            opts.splitting.eval_matrix_functions(eqn, opts, oper, xj(j));
       if eqn.type == 'T'
           LB{j} = eqn.B;
       elseif eqn.type == 'N'
           LB{j} = eqn.C';
       end
    end
    LB = cell2mat(LB);
    
    if eqn.type == 'T'
        DB = kron(diag(wj), eqn.Rinv*speye(size(LB,2) / length(xj)));
    elseif eqn.type == 'N'
        DB = kron(diag(wj), speye(size(LB,2) / length(xj)));
    end
%     if eqn.type == 'T'
%         L = cell2mat(arrayfun(eqn.B, xj', 'UniformOutput', false));
%         D = kron(diag(wj), eqn.Rinv*speye(size(L,2) / length(xj)));
%     elseif eqn.type == 'N'
%         L = cell2mat(arrayfun(eqn.C', xj', 'UniformOutput', false));
%         D = kron(diag(wj), speye(size(L,2) / length(xj)));
%     end
            
    [LB, DB] = mess_column_compression(LB, 'N', DB, ...
        opts.splitting.trunc_tol, opts.splitting.trunc_info);

	D = (I + D0*L0'*(LB*(DB*(LB'*L0))))\D0;
end

end