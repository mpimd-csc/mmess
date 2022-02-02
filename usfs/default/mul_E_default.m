function C=mul_E_default(eqn, opts,opE,B,opB)

% function C=mul_E_default(eqn, opts,opE,B,opB)
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
% This function uses another 'default' usfs function
% (size_default(eqn,opts)) to obtain the number of rows of matrix A_
% in structure eqn, that should be equal to the number of rows of
% the matrix E_.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input parameters
if (not(ischar(opE)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opE or opB is not a char');
end

opE = upper(opE); opB = upper(opB);
if not((opE=='N' || opE=='T'))
    error('MESS:error_arguments', 'opE is not ''N'' or ''T''');
end

if not((opB=='N' || opB=='T'))
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments', 'B has to ba a matrix');
end

%% check data in eqn structure
if not(isfield(eqn,'E_'))
    error('MESS:error_arguments', 'field eqn.E_ is not defined');
end

rowE = size_default(eqn, opts);
colE = rowE;

%% perform multiplication
switch opE

    case 'N'
        switch opB

            % implement operation E_ * B
            case 'N'
                if not(colE == size(B,1))
                    error('MESS:error_arguments', ...
                          ['number of columns of E_ differs with number ' ...
                           'of rows of B']);
                end
                C = eqn.E_ * B;

            % implement operation E_ * B'
            case 'T'
                if not(colE == size(B,2))
                    error('MESS:error_arguments', ...
                          ['number of columns of E_ differs with number ' ...
                           'of columns of B']);
                end
                C = eqn.E_ * B';
        end

    case 'T'
        switch opB

            % implement operation E_' * B
            case 'N'
                if not(rowE == size(B,1))
                    error('MESS:error_arguments', ...
                          ['number of rows of E_ differs with number ' ...
                           'of rows of B']);
                end
                C = eqn.E_' * B;

            % implement operation E_' * B'
            case 'T'
                if not(rowE == size(B,2))
                    error('MESS:error_arguments', ...
                          ['number of rows of E_ differs with number ' ...
                           'of columns of B']);
                end
                C = eqn.E_' * B';
        end

end

end
