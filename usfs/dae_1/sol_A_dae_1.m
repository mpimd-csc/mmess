function X = sol_A_dae_1(eqn, opts, opA, B, opB)%#ok<INUSL>
%  function sol_A solves solves opA(A_)*X = opB(B)
%
% Input:
%  eqn       structure with field A_
%  opts      struct contains parameters for the algorithm
%  opA       character specifies the form of \opA(A_)
%                  opA = 'N' solves A_*X=opB(B)
%                  opA = 'T' solves A_^T*X=opB(B)
%
%  B         p-x-q matrix
%
%  opB       character specifies the form of opB(B)
%                  opB = 'N' solves A_*X=B
%                  opB = 'T' solves A_*X=B^T
%
% Output:
%  X       matrix fulfills equation opA(A_)X = opB(B)
%
% uses size_dae_1

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%



%% check input Parameters
if (not(ischar(opA)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opA or opB is not a char');
end

opA = upper(opA); opB = upper(opB);
if(not((opA == 'N' || opA == 'T')))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(not(isfield(eqn, 'A_'))) || not(isnumeric(eqn.A_))
    error('MESS:error_arguments', 'field eqn.A_ is not defined');
end

if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or corrupted st field detected in equation structure.');
end

n = size(eqn.A_,1);
st = eqn.st;
[rowB,colB] = size(B);

if(opB == 'N')
    if(n > rowB)
        B = [B;
            zeros(n - st, colB)];
    elseif n < rowB
        error('MESS:error_arguments', 'B has more rows than A');
    end
else
    if(n > colB)
        B = [B, zeros(rowB, n - st)];
    elseif n < colB
        error('MESS:error_arguments', 'B has more columns than A');
    end
end

%% solve
switch opA

    case 'N'
        switch opB

            %implement solve A_*X=B
            case 'N'
                if(n ~= size(B, 1))
                    error('MESS:error_arguments', ...
                          'number of rows of A_ differs with rows of B');
                end
                X = eqn.A_ \ B;

            %implement solve A_*X=B'
            case 'T'
                if(n ~= size(B, 2))
                    error('MESS:error_arguments', ...
                          'number of rows of A_ differs with cols of B');
                end
                X = eqn.A_ \ B';
        end

    case 'T'
        switch opB

            %implement solve A_'*X=B
            case 'N'
                if(n ~= size(B, 1))
                    error('MESS:error_arguments', ...
                          'number of cols of A_ differs with rows of B');
                end
                X = eqn.A_' \ B;

            %implement solve A_'*X=B'
            case 'T'
                if(n ~= size(B, 2))
                    error('MESS:error_arguments', ...
                          'number of cols of A_ differs with cols of B');
                end
                X = eqn.A_' \ B';
        end

end
X = X(1 : st, :);
end
