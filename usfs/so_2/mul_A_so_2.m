function C=mul_A_so_2(eqn, opts,opA,B,opB)%#ok<INUSL>
% function C=mul_A_so_2(eqn, opts,opA,B,opB)
%
% Call help mess_usfs_so_2 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns C = A*B, where matrix A given by structure eqn and input matrix B could be transposed.
% Matrix A is assumed to be quadratic and has a size of 2* size(K).
%
%   Inputs:
%
%   eqn     structure containing  data for matrix A (fields 'K_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' performs A*opB(B)
%           opA = 'T' performs A'*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs opA(A)*B
%           opB = 'T' performs opA(A)*B'
%
%   Output:
%
%       |-K  0|
%   C = | 0  M| *opB(B)
%
% This function does not use other so3 functions.
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
if(not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_)) || not(isfield(eqn,'M_')) ...
        || not(isnumeric(eqn.M_)))
    error('MESS:error_arguments',...
        'A consists of K and M, field eqn.K_ or eqn.M_ is not defined');
end


[rowK, colK] = size(eqn.K_);

colA = 2*colK;

%% perform multiplication
switch opB

    % implement operation A*B = C
    case 'N'
        if(colA ~= size(B,1))
              error('MESS:error_arguments','number of columns of A differs with number of rows of B');
        end

        C = [-eqn.K_*B(1:rowK,:) ;...
              eqn.M_*B(rowK+1:end,:)];

    % implement operation A*B' = C
    case 'T'
        if(colA ~= size(B,2))
            error('MESS:error_arguments','number of columns of A differs with number of columns of B');
        end

        C = [-eqn.K_*B(:,1:colK)';...
              eqn.M_*B(:,colK+1:end)'];
end

end
