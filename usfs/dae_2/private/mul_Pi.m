function C = mul_Pi(eqn, opts, type, opP, B, opB)
% MUL_Pi multiplies with the hidden manifold projection matrix or its
% transpose. Note that the multiplication is actually implemented as the
% solution of a saddle point structured linear system.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isfield(eqn, 'P_'))
    mess_err(opts, 'equation_data', ...
             'Could not find required field ''P_''. Did you not run mul_Pi_pre?');
end

n      = size(eqn.P_, 1);
one    = 1:eqn.manifold_dim;
n_two  = n - eqn.manifold_dim;

%% build extended RHS
switch opB
    case 'N'
        H = B;
    case 'T'
        H = B';
    otherwise
        mess_err(opts, 'input_data', ...
                 'opB must be either ''N'' or ''T''.');
end

%% solve augmented saddle point system
if type == 'r'
    if opP == 'N'
        H = [eqn.E_(one, one) * H; zeros(n_two, size(H, 2))];
        C = eqn.P_ \ H;
        C = C(one, :);
    elseif opP == 'T'
        H = [H; zeros(n_two, size(H, 2))];
        C = eqn.P_' \ H;
        C = eqn.E_(one, one)' * C(one, :);
    else
        mess_err(opts, 'input_data', ...
                 'opP must be either ''N'' or ''T''.');

    end
elseif type == 'l'
    if opP == 'N'
        H = [H; zeros(n_two, size(H, 2))];
        C = eqn.P_ \ H;
        C = eqn.E_(one, one) * C(one, :);
    elseif opP == 'T'
        H = [eqn.E_(one, one)' * H; zeros(n_two, size(H, 2))];
        C = eqn.P_' \ H;
        C = C(one, :);
    else
        mess_err(opts, 'input_data', ...
                 'opP must be either ''N'' or ''T''.');
    end
else
    mess_err(opts, 'input_data', ...
             'type must be either ''N'' or ''T''.');
end

end
