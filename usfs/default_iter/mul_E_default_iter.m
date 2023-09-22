function C = mul_E_default_iter(eqn, opts, opE, B, opB)

% function C=mul_E_default_iter(eqn, opts,opE,B,opB)
%
% This function returns C = E_*B, where matrix E_ given by structure
% eqn and input matrix B could be transposed. Matrix E_ is assumed
% to be quadratic and has the same size as A_ in structure eqn.
%
%   Inputs:
%
%   eqn     structure containing field 'E_'
%   opts    structure containing parameters for the algorithm
%   opE     character specifying the shape of E_
%           opE = 'N' performs E_*opB(B)
%           opE = 'T' performs E_'*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs opE(E_)*B
%           opB = 'T' performs opE(E_)*B'
%
%   Output:
%
%   C = opE(E_)*opB(B)
%
% This function uses another default function
% size_default_iter(eqn,opts) to obtain the number of rows of matrix A_
% in structure eqn, that should be equal to the number of rows of
% the matrix E_.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check input parameters
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
    mess_err(opts, 'error_arguments', 'B has to be a matrix');
end

%% Check data in eqn structure
if not(isfield(eqn, 'E_'))
    mess_err(opts, 'error_arguments', 'field eqn.E_ is not defined');
end

rowE = size_default_iter(eqn, opts);
colE = rowE;

%% Perform multiplication
switch opE

    case 'N'
        switch opB

            % Implement multiplication E_*B
            case 'N'
                if not(colE == size(B, 1))
                    mess_err(opts, 'error_arguments', ...
                             ['number of columns of E_ differs from number ' ...
                              'of rows of B']);
                end
                C = eqn.E_ * B;

                % Implement multiplication E_*B'
            case 'T'
                if not(colE == size(B, 2))
                    mess_err(opts, 'error_arguments', ...
                             ['number of columns of E_ differs from number ' ...
                              'of columns of B']);
                end
                C = eqn.E_ * B';
        end

    case 'T'
        switch opB

            % Implement multiplication E_'*B
            case 'N'
                if not(rowE == size(B, 1))
                    mess_err(opts, 'error_arguments', ...
                             ['number of rows of E_ differs from number ' ...
                              'of rows of B']);
                end
                C = eqn.E_' * B;

                % Implement multiplication E_'*B'
            case 'T'
                if not(rowE == size(B, 2))
                    mess_err(opts, 'error_arguments', ...
                             ['number of rows of E_ differs from number ' ...
                              'of columns of B']);
                end
                C = eqn.E_' * B';
        end

end

end
