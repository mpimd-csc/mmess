function C = mul_Pi(eqn,opP,B, opB)
% MUL_Pi multiplies with the hidden manifold projection matrix or its
% transpose. Note that the multiplication is actually implemented as the
% solution of a saddle point structured linear system.
%

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


if not(isfield(eqn,'P_'))
    error('MESS:equation_data',['Could not find required field ''P_''.',...
        'Did you not run mul_Pi_pre? ']);
end

n = size(eqn.P_,1);
st = eqn.st;

if opP == 'N'
    if opB == 'N'
        C = eqn.P_ \ [B; zeros(n -st,size(B,2))];
        C = eqn.E_(1:eqn.st,1:eqn.st) * C(1:eqn.st,:);
    elseif opB== 'T'
        C = eqn.P_ \ [B, zeros(size(B,1),n-st)]';
        C = eqn.E_(1:eqn.st,1:eqn.st) * C(1:eqn.st,:);
    else
        error('MESS:input_data','opB must be either ''N'' or ''T''.');
    end
elseif opP == 'T'
    if opB == 'N'
        H = B;
    elseif opB== 'T'
        H = B';
    else
        error('MESS:input_data','opB must be either ''N'' or ''T''.');
    end
    H = [eqn.E_(1:eqn.st,1:eqn.st)' * H; zeros(n-st,size(H,2))];
    C = eqn.P_' \ H;
    C = C(1:eqn.st,:);

else
    error('MESS:input_data','opP must be either ''N'' or ''T''.');
end

end

