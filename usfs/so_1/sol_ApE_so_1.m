function X=sol_ApE_so_1(eqn, opts,opA,p,opE,C,opC)%#ok<INUSL>
% function X=sol_ApE_so_1(eqn, opts,opA,p,opE,C,opC)
%
% Call help mess_usfs_so_1 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns X =(A+p*E)\C, where matrices A and E given
% by structure eqn and input matrix C could be transposed. Matrices
% A and E are assumed to be quadratic.
%
%   Inputs:
%
%   eqn     structure containing data for matrices A
%           (fields 'K_' and 'E_') and E (fields 'K_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' solves (A + p* opE(E))*X = opC(C)
%           opA = 'T' solves (A' + p* opE(E))*X = opC(C)
%   p       scalar value
%   opE     character specifying the shape of E
%           opE = 'N' solves (opA(A) + p* E)*X = opC(C)
%           opE = 'T' solves (opA(A) + p* E')*X = opC(C)
%   C       n-x-p matrix
%   opC     character specifies the shape of C
%           opC = 'N' solves (opA(A) + p* opE(E))*X = C
%           opC = 'T' solves (opA(A) + p* opE(E))*X = C'
%
%   Output:
%
%                                       (| 0 -K|     |-K  0|)
%   X       matrix fullfilling equation (|-K  E| + p*| 0  M|)*X = opC(C)
%
%   This function does not use other so1 functions.
%
% ATTENTION: opA and opE are not used since matrices A and E are symmetric

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input parameters
if (not(ischar(opA)) || not(ischar(opE)) || not(ischar(opC)))
    error('MESS:error_arguments', 'opA, opE or opC is not a char');
end

opA = upper(opA); opE = upper(opE); opC = upper(opC);

if(not((opA=='N' || opA=='T')))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(not((opE=='N' || opE=='T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not((opC=='N' || opC=='T')))
    error('MESS:error_arguments','opC is not ''N'' or ''T''');
end

if(not(isnumeric(p)))
    error('MESS:error_arguments','p is not numeric');
end

if (not(isnumeric(C))) || (not(ismatrix(C)))
    error('MESS:error_arguments','C has to ba a matrix');
end

[rowK, colK] = size(eqn.K_);

colA = 2*colK;
one = 1:rowK;
two = (rowK + 1) : colA;

switch opC

    case 'N'    % implement solve (A+p*E)*X=C

        if not(colA == size(C,1))
            error('MESS:error_arguments',['number of rows of A ' ...
                'differs with number of rows of C']);
        end

        X2 = (p * (p * eqn.M_ - eqn.E_) + eqn.K_) \ ...
             (p * C(two,:) - C(one,:));
        X1 = (-p*eqn.K_) \ (C(one,:) + eqn.K_ * X2);
        X = [ X1; X2];

    case 'T'   % implement solve (A+p*E)*X=C'

        if not(colA == size(C,2))
            error('MESS:error_arguments',['number of rows of A ' ...
                'differs with number of columns of C']);
        end

        X2 = (p * (p * eqn.M_ - eqn.E_) + eqn.K_) \ ...
             (p * C(:, two)' - C(:, one)');
        X1 = -(p * eqn.K_) \ (C(:, one)' + eqn.K_ * X2);
        X = [X1; X2];
end
end

