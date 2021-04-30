function X=sol_ApE_so_2(eqn, opts,opA,p,opE,B,opB)%#ok<INUSL>
% function X=sol_ApE_so_2(eqn, opts,opA,p,opE,C,opC)
%
% Call help mess_usfs_so_2 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns X =(A+p*E)\C, where matrices A and E given
% by structure eqn and input matrix C could be transposed.
% Matrices A and E are assumed to be quadratic.
%
%   Inputs:
%
%   eqn     structure containing data for matrices
%           A (fields 'K_' and 'M_') and E (fields 'E_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' solves (A + p* opE(E))*X = opC(C)
%           opA = 'T' solves (A' + p* opE(E))*X = opC(C)
%   p       scalar value
%   opE     character specifying the shape of E
%           opE = 'N' solves (opA(A) + p* E)*X = opC(C)
%           opE = 'T' solves (opA(A) + p* E')*X = opC(C)
%   B       n-x-p matrix
%   opB     character specifies the shape of B
%           opC = 'N' solves (opA(A) + p* opE(E))*X = B
%           opC = 'T' solves (opA(A) + p* opE(E))*X = B'
%
%   Output:
%
%                                       (|-K  0|     |E  M|)
%   X       matrix fullfilling equation (| 0  M| + p*|M  0|)*X = opB(B)
%
%   This function does not use other so3 functions.
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
if (not(ischar(opA)) || not(ischar(opE)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opA, opE or opC is not a char');
end

opA = upper(opA); opE = upper(opE); opB = upper(opB);

if(not((opA=='N' || opA=='T')))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(not((opE=='N' || opE=='T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not((opB=='N' || opB=='T')))
    error('MESS:error_arguments','opC is not ''N'' or ''T''');
end

if(not(isnumeric(p)))
    error('MESS:error_arguments','p is not numeric');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
if(eqn.haveE ==1)
    if(not(isfield(eqn,'M_')) || not(isnumeric(eqn.M_)) || not(isfield(eqn,'E_')) || ...
            not(isnumeric(eqn.E_)) || not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_)))
        error('MESS:error_arguments',...
            'field eqn.M_, eqn.E_ or eqn.K_ is not defined or corrupted');
    end
else
    error('MESS:error_arguments',['eqn.haveE has to be 1 because of ' ...
                        'the structure of E']);
end

rowK = size(eqn.K_, 1);
rowA = 2*rowK;

%% perform solve operations
switch opB
    % implement solve (A+p*E)*X = B
    case 'N'
        if (rowA ~= size(B,1))
            error('MESS:error_arguments',['number of rows of A differs ' ...
                                'with number of rows of B']);
        end
        temp = p*eqn.E_ - eqn.K_;
        X1 = (p^2*eqn.M_ - temp)\(p*B(rowK+1:end,:) - B(1:rowK,:));
        X2 = (p*eqn.M_)\(B(1:rowK,:) - temp*X1);
        X  = [X1;X2];

    % implement solve (A+p*E)*X = B'
    case 'T'
        if (rowA ~= size(B,2))
            error('MESS:error_arguments',['number of rows of A differs ' ...
                                'with number of columns of B']);
        end
        temp = p*eqn.E_ -eqn.K_;
        X1 = (p^2*eqn.M_ - temp)\(p*B(:,rowK+1:end)' - B(:,1:rowK)');
        X2 = (p*eqn.M_)\(B(:,1:rowK)' - temp*X1);
        X  = [X1;X2];
end
