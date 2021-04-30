function [rw,  Hp, Hm, Vp, Vm, eqn, opts, oper] = get_ritz_vals_so_2(eqn, opts, oper, U, W, p_old)
% [rw,  Hp, Hm, Vp, Vm, eqn, opts, oper] = get_ritz_vals_so_2(eqn,opts,oper)
%
% Call help mess_usfs_so_2 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns suitable Ritz values, Hessenberg matrices
% and matrices consisting of basis vectors corresponding to A and
% A^{-1}.
%
%   Input:
%
%   eqn      data structure
%   opts     structure containing parameters for the algorithm
%   oper
%
%   Output:
%
%   rw       vector of Ritz values
%   Hp       Hessenberg matrix corresponding to A
%   Hm       Hessenberg matrix corresponding to A^{-1}
%   Vp       matrix consisting of basis vectors corresponding to A
%   Vm       matrix consisting of basis vectors corresponding to A^{-1}
%
% This function does not use other so3 functions.
% This function uses another help function
% mess_get_ritz_vals(eqn,opts,oper), which returns all Ritz values
% and Hessenberg matrices and matrices consisting of basis vectors
% corresponding to A and A^{-1}.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


if isfield(opts.shifts, 'method') && ...
        strcmp(opts.shifts.method, 'projection')
    if isempty(W)
        % first shifts are computed with U = eqn.G and W = A * eqn.G
        W = oper.mul_A(eqn, opts, eqn.type, U, 'N');
    end
    rw = mess_projection_shifts(eqn, opts, oper, U, ...
        W, p_old);
else
    n = oper.size(eqn, opts);
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
if isfield(opts.shifts,'truncate') && isnumeric(opts.shifts.truncate)
    rw = rw(abs(rw) < opts.shifts.truncate);
    rw = rw(abs(rw) > 1/opts.shifts.truncate);
end
end
