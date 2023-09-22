function C = mul_E_dae_2(eqn, opts, opE, B, opB)

%% function mul_E performs operation C = opE(E_)*opB(B)
%
% Input:
%   eqn     structure contains field E_
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' performs E_*opB(B)
%           opE = 'T' performs E_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opE(E_)*B
%           opB = 'T' performs opE(E_)*B'
%
% Output:
% C = opE(E_)*opB(B)
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
if (not(isfield(eqn, 'E_'))) || not(isnumeric(eqn.E_))
    mess_err(opts, 'error_arguments', 'field eqn.E_ is not defined');
end

if (not(isfield(eqn, 'M_'))) || not(isnumeric(eqn.M_))
    mess_err(opts, 'error_arguments', ...
             'field eqn.M_ is not defined. Did you forget to run mul_E_pre?');
end
if not(isfield(eqn, 'manifold_dim'))    || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'error_arguments', ...
             ['Missing or Corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end

n = size(eqn.E_, 1);
one = 1:eqn.manifold_dim;

switch opB
    case 'N'
        dim = size(B, 1);
    case 'T'
        dim = size(B, 2);
end

if not(dim == eqn.manifold_dim) && not(dim == n)
    mess_err(opts, 'error_arguments', 'size of B does not match data in E');
end

%% perform multiplication
if dim == n
    switch opE

        case 'N'

            switch opB

                case 'N'
                    % implement operation A_*B
                    C = eqn.M_ * B;

                case 'T'
                    % implement operation A_*B'
                    C = eqn.M_ * B';

            end

        case 'T'

            switch opB

                case 'N'
                    % implement operation A_'*B
                    C = eqn.M_' * B;

                case 'T'
                    % implement operation A_'*B'
                    C = eqn.M_' * B';

            end

    end

else

    switch opE
        case 'N'
            V = eqn.E_(one, one) * mul_Pi(eqn, opts, 'r', 'N',  B, opB);
            C = mul_Pi(eqn, opts, 'l', 'N', V, 'N');
        case 'T'
            V = eqn.E_(one, one)' * mul_Pi(eqn, opts, 'l', 'T',  B, opB);
            C = mul_Pi(eqn, opts, 'r', 'T', V, 'N');
    end
end
end
