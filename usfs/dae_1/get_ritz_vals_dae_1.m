function [rw, Hp, Hm, Vp, Vm] = get_ritz_vals_dae_1(eqn, opts, oper, U, W, p_old)
%  This function ensures that W is not empty if the shift
%  method is projection. Otherwise, it checks opts.shifts.b0. It
%  should be a vector of the same size as eqn.A but oper.size gives eqn.st.

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
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
