function X=sol_A_so_1(eqn, opts,opA,B,opB)%#ok<INUSL>
% function X=sol_A_so_1(eqn, opts,opA,B,opB)
%
% Call help mess_usfs_so_1 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns X = A\B, where matrix A given by structure eqn and input matrix B could be transposed.
% Matrix A is assumed to be quadratic and has a size of 2*size(K).
%
%   Inputs:
%
%   eqn     structure containing data for matrix A (fields 'E_' and 'K_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' solves A *X = opB(B)
%           opA = 'T' sovles A'*X = opB(B)
%   B       p-x-q matrix
%   opB     character specifying the shape of B
%           opB = 'N' solves opA(A)*X = B
%           opB = 'T' solves opA(A)*X = B'
%
%   Output:
%                                       | 0 -K|
%   X       matrix fullfilling equation |-K -E| *X= opB(B)
%
%   This function does not use other so1 functions.
%
% ATTENTION: opA is not used since matrix A is symmetric

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input parameters
if (not(ischar(opA)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opA or opB is not a char');
end

opA = upper(opA); opB = upper(opB);
if(not((opA=='N' || opA=='T')))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(not((opB=='N' || opB=='T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_)) || not(isfield(eqn,'E_'))) || not(isnumeric(eqn.E_))
    error('MESS:error_arguments',...
        'A consists of K and D, field eqn.K_ or eqn.E_ is not defined or corrupted');
end

[rowK, colK] = size(eqn.K_);
rowA = 2*rowK;



%% perform solve operations
switch opB

    %implement solve A*X=B
    case 'N'

        if(rowA~=size(B,1))
            error('MESS:error_arguments','number of rows of A differs with number of rows of B');
        end

        %%K_ hat vollen Rang
        X2 = eqn.K_\B(1:rowK,:);
        X1 = eqn.K_\(B(rowK+1:end,:)-eqn.E_*X2);
        X= [-X1;-X2];

    %implement solve A*X=B'
    case 'T'

        if(rowA~=size(B,2))
            error('MESS:error_arguments','number of rows of A differs with number of columns of B');
        end

        X2 = eqn.K_\B(:,1:colK)';
        X1 = eqn.K_\(B(:,colK+1:end)'-eqn.E_*X2);
        X= [-X1;-X2];
end


end

