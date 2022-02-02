function X=sol_A_default(eqn, opts,opA,B,opB)
% function X=sol_A_default(eqn, opts,opA,B,opB)
%
% This function returns X = A_\B, where matrix A_ given by
% structure eqn and input matrix B could be transposed. Matrix A_
% is assumed to be quadratic.
%
%    Inputs:
%
%    eqn       structure containing field 'A_'
%    opts      structure containing parameters for the algorithm
%    opA       character specifying the shape of A_
%                  opA = 'N' solves  A_*X = opB(B)
%                  opA = 'T' solves  A_'*X = opB(B)
%    B         p-x-q matrix
%    opB       character specifying the shape of B
%                  opB = 'N' solves  opA(A_)*X = B
%                  opB = 'T' solves  opA(A_)*X = B'
%
%    Output:
%
%    X         matrix fulfilling equation  opA(A_)*X = opB(B)
%
% This function uses another default function size_default(eqn,
% opts) to obtain the number of rows of matrix A_ in structure
% eqn.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%



%% check input parameters
if not(ischar(opA)) || not(ischar(opB))
    error('MESS:error_arguments', 'opA or opB is not a char');
end

opA = upper(opA);
opB = upper(opB);
if not( opA=='N' || opA=='T' )
    error('MESS:error_arguments', 'opA is not ''N'' or ''T''');
end

if not( opB=='N' || opB=='T' )
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end

if not( isnumeric(B)) ||  not(ismatrix(B))
    error('MESS:error_arguments', 'B has to ba a matrix');
end

%% check data in eqn structure
if not(isfield(eqn,'A_'))
    error('MESS:error_arguments', 'field eqn.A_ is not defined');
end

rowA = size_default(eqn, opts);
colA = rowA;

%% perform solve operations
switch opA

  case 'N'
    switch opB

      %implement solve A_ * X = B
      case 'N'
        if not(rowA == size(B,1))
            error('MESS:error_arguments', ...
                  'number of rows of A_ differs with number of rows of B');
        end
        X = eqn.A_ \ B;

      %implement solve A_ * X = B'
      case 'T'
        if not(rowA == size(B,2))
            error('MESS:error_arguments', ...
                  'number of rows of A_ differs with number of columns of B');
        end
        X = eqn.A_ \ B';
    end

  case 'T'
    switch opB

      %implement solve A_' * X = B
      case 'N'
        if not(colA == size(B,1))
            error('MESS:error_arguments', ...
                  'number of columns of A_ differs with number of rows of B');
        end
        X = eqn.A_' \ B;

        %implement solve A_' * X = B'
      case 'T'
        if not(colA == size(B,2))
            error('MESS:error_arguments', ['number of columns of A_ ' ...
                                'differs with number of columns of B']);
        end
        X = eqn.A_' \ B';
    end

end

end

