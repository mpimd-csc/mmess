function X=sol_E_so_1(eqn, opts,opE,B,opB)%#ok<INUSL>
% function X=sol_E_so_1(eqn, opts,opE,B,opB)
%
% Call help mess_usfs_so_1 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns X= E\B, where matrix E given by structure eqn and input matrix B could be transposed.
%
%   Inputs:
%   eqn     structure containing data for matrix E (fields 'M_' and 'K_')
%   opts    structure containing parameters for the algorithm
%   opE     character specifying the shape of E
%           opE = 'N' solves E *X = opB(B)
%           opE = 'T' sovles E'*X = opB(B)
%   B       p-x-q matrix
%   opB     character specifying the shape of B
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
%   Output:
%
%                                       |-K  0|
%   X       matrix fullfilling equation | 0  M| *X = opB(B)
%
%   This function does not use other so1 functions.
%
% ATTENTION: opE is not used since matrix E is symmetric

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%



%% check input parameters
if (not(ischar(opE)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opE or opB is not a char');
end
opE = upper(opE); opB = upper(opB);
if(not((opE=='N' || opE=='T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not((opB=='N' || opB=='T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_)) || not(isfield(eqn,'M_')) ...
        || not(isnumeric(eqn.M_)))
    error('MESS:error_arguments',...
        'E consists of K and M, field eqn.K_ or eqn.M_ is not defined');
end


[rowK, colK] = size(eqn.K_);
rowE = 2*rowK;


%% perform solve operations
switch opB

    %implement solve E*X=B
    case 'N'
        if(rowE~=size(B,1))
            error('MESS:error_arguments','number of rows of E differs with number of rows of B');
        end
        X= [-eqn.K_\B(1:rowK,:);eqn.M_\B(rowK+1:end,:)];

    %implement solve E*X=B'
    case 'T'
        if(rowE~=size(B,2))
            error('MESS:error_arguments','number of rows of E differs with number of columns of B');
        end
        X= [-eqn.K_\B(:,1:colK)';eqn.M_\B(:,colK+1:end)'];
end

end

