function X = sol_E_default(eqn, opts, opE, B, opB)
% function X=sol_E_default(eqn, opts,opE,B,opB)
%
% This function returns X = E_\B, where matrix E_ given by
% structure eqn and input matrix B could be transposed. Matrix E_
% is assumed to be quadratic and has the same size as A_ in
% structure eqn.
%
%   Inputs:
%
%   eqn       structure containing field 'E_'
%   opts      structure containing parameters for the algorithm
%   opE       character specifying the shape of E_
%                  opE = 'N' solves E_*X = opB(B)
%                  opE = 'T' solves E_'*X = opB(B)
%   B         p-x-q matrix
%   opB       character specifying the shape of B
%                  opB = 'N' solves  opE(E_)*X = B
%                  opB = 'T' solves  opE(E_)*X = B'
%
%   Output:
%
%   X         matrix fulfilling equation  opE(E_)*X = opB(B)
%
% This function uses another default function size_default(eqn,
% opts) to obtain the number of rows of matrix A_ in structure eqn,
% that should be equal to the number of rows of matrix E_.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input parameters
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

if not(isnumeric(B)) || not(ismatrix(B))
    mess_err(opts, 'error_arguments', 'B has to ba a matrix');
end

%% check data in eqn structure
if not(isfield(eqn, 'E_'))
    mess_err(opts, 'error_arguments', 'field eqn.E_ is not defined');
end

rowE = size_default(eqn, opts);
colE = rowE; % we only consider square systems

%% perform solve operations
switch opE

    case 'N'
        switch opB

            % implement solve E_ * X = B
            case 'N'
                if not(rowE == size(B, 1))
                    mess_err(opts, 'error_arguments', ...
                             ['number of rows of E_ differs with number ' ...
                              'of rows of B']);
                end
                X = eqn.E_ \ B;

                % implement solve E_ * X = B'
            case 'T'
                if not(rowE == size(B, 2))
                    mess_err(opts, 'error_arguments', ...
                             ['number of rows of E_ differs with number ' ...
                              'of columns of B']);
                end
                X = eqn.E_ \ B';
        end

    case 'T'
        switch opB

            % implement solve E_' * X = B
            case 'N'
                if not(colE == size(B, 1))
                    mess_err(opts, 'error_arguments', ...
                             ['number of columns of E_ differs with number ' ...
                              'of rows of B']);
                end
                X = eqn.E_' \ B;

                % implement solve E_' * X = B'
            case 'T'
                if not(colE == size(B, 2))
                    mess_err(opts, 'error_arguments', ['number of columns ' ...
                                                       'of E_ differs with number ' ...
                                                       'of columns of B']);
                end
                X = eqn.E_' \ B';
        end

end

end
