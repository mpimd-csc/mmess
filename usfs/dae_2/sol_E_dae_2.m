function X = sol_E_dae_2(eqn, opts, opE, B, opB)
%% function sol_E solves opE(M_)*X = opB(B) resp. performs X=opE(M_)\opB(B)
% sol_E_pre should be called before to construct
% M_ = [ E1 -J';
%      [ J   0 ]
% from
% A = [ A1 -J';
%       J   0]
% E = [ E1 0;
%       0  0]
%
% Input:
%   eqn     structure contains data for M_
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' solves E *X = opB(B)
%           opE = 'T' solves E'*X = opB(B)
%
%   B       p-x-q matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
% Output
%
%   X       matrix fulfills equation opE(E)*X = opB(B)
%
%   uses no other dae_2 function

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input Parameters
if not(ischar(opE)) || not(ischar(opB))
    mess_err(opts, 'error_arguments', 'opE or opB is not a char');
end

opE = upper(opE);
opB = upper(opB);
if not(opE == 'N' || opE == 'T')
    mess_err(opts, 'error_arguments', 'opE is not ''N'' or ''T''');
end

if not(opB == 'N' || opB == 'T')
    mess_err(opts, 'error_arguments', 'opB is not ''N'' or ''T''');
end
if (not(isnumeric(B))) || (not(ismatrix(B)))
    mess_err(opts, 'error_arguments', 'B has to ba a matrix');
end

%% check data in eqn structure
if not(isfield(eqn, 'M_'))
    mess_err(opts, 'error_arguments', ['field eqn.M_ is not defined. Did ' ...
                                       'you forget to run sol_E_pre?']);
end
if not(isfield(eqn, 'manifold_dim'))    || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equation_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end

n = size(eqn.M_, 1);

[rowB, colB] = size(B);

if opB == 'N'
    if rowB == eqn.manifold_dim
        B = [B; zeros(n - rowB, colB)];
    elseif not(rowB == n)
        mess_err(opts, 'error_arguments', 'size of B does not match data in E');
    end
else
    if colB == eqn.manifold_dim
        B = [B, zeros(rowB, n - colB)];
    elseif not(colB == n)
        mess_err(opts, 'error_arguments', 'size of B does not match data in E');
    end
end

%% solve
switch opE

    case 'N'
        switch opB

            % implement solve M_*X=B
            case 'N'

                X = eqn.M_ \ B;

                % implement solve M_*X=B'
            case 'T'

                X = eqn.M_ \ B';
        end

    case 'T'
        switch opB

            % implement solve M_'*X=B
            case 'N'

                X = eqn.M_' \ B;

                % implement solve M_'*X=B'
            case 'T'

                X = eqn.M_' \ B';
        end

end

end
