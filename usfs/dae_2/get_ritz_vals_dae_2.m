function [rw, Hp, Hm, Vp, Vm] = get_ritz_vals_dae_2(eqn, opts, oper, U, W, p_old)
%  This function is an exact copy of the Penzl heuristic part in mess_para.
%  the only difference is that B or C and K are filled up by trailing zero
%  blocks to allow for the computation of the Ritz values with respect to
%  the full size block matrices instead of the restriction to the
%  (1,1)-block for which the ADI is formulated.
%
%  MMESS (Jens Saak, October 2013)

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
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
% returns order of A or states of A, A is supposed to be square
n = size(eqn.A_,1);

%%
% here we add the trailing zero blocks. Note that we are not passing the
% eqn structure back as an output, so this change is not visible in
% anything above this routine and will only be passed on to the function
% handles used in here.
if isfield(eqn, 'U') && not(isempty(eqn.U))
    eqn.U = [eqn.U; zeros(n - size(eqn.U, 1), size(eqn.U, 2))];
end
if isfield(eqn,'V') && not(isempty(eqn.V))
    eqn.V = [eqn.V; zeros(n - size(eqn.V,1), size(eqn.V,2))];
end
if isfield(opts.shifts, 'method') && ...
        strcmp(opts.shifts.method, 'projection')
    U = [U; zeros(n - size(U, 1), size(U, 2))];
    if isempty(W)
        % first shifts are computed with U = eqn.G and W = A * eqn.G
        W = oper.mul_A(eqn, opts, eqn.type, U, 'N');
    else
        W = [W; zeros(n - size(W, 1), size(W, 2))];
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
