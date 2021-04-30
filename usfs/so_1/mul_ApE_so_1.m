function C=mul_ApE_so_1(eqn, opts,opA,p,opE,B,opB)%#ok<INUSL>
% function C=mul_ApE_so_1(eqn, opts,opA,p,opE,B,opB)
%
% Call help mess_usfs_so_1 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns C = (A+p*E)*B, where matrices A and E given
% by structure eqn and input matrix B could be transposed.
% Matrix A is assumed to be quadratic and has a size of 2*size(K).
%
%   Inputs:
%
%   eqn     structure containing data for matrices A
%           (fields 'K_' and 'E_') and E (fields 'K_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' performs (A + p*opE(E))*opB(B)
%           opA = 'T' performs (A' + p*opE(E))*opB(B)
%   p       scalar value
%   opE     character specifying the shape of E
%           opA = 'N' performs (opA(A) + p*E)*opB(B)
%           opA = 'T' performs (opA(A) + p*E')*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs (opA(A) + p*opE(E))*B
%           opB = 'T' performs (opA(A) + p*opE(E))*B'
%
%   Output:
%
%       (| 0 -K|      |-K  0|)
%   C = (|-K -E| + p* | 0  M|)*opB(B)
%
% This function does not use other so1 functions.
%
% ATTENTION: opA, opE are not used since matrices A and E are symmetric.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input parameters
if (not(ischar(opA)) || not(ischar(opB)) || not(ischar(opE)))
    error('MESS:error_arguments', 'opA, opB or opE is not a char');
end

opA = upper(opA); opB = upper(opB); opE = upper(opE);

if(not((opA=='N' || opA=='T')))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(not((opB=='N' || opB=='T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if(not((opE=='N' || opE=='T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not(isnumeric(p)))
    error('MESS:error_arguments','p is not numeric');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure

[rowK, colK] = size(eqn.K_);
colA = 2*colK;
one = 1:rowK;
two = (rowK + 1) : colA;

switch opB

      %implement multiplication (A+p*E)*B=C
      case 'N'
        if(colA~=size(B,1))
            error('MESS:error_arguments',['number of columns of A ' ...
                                'differs with number of rows of B']);
        end
        C = -[(p*(eqn.K_*B(one, :)) + eqn.K_*B(two, :));
              (eqn.K_*B(one, :)) - p*(eqn.M_*B(two, :)) + eqn.E_*B(two, :)];

      %implement multiplication (A+p*E)*B'=C
      case 'T'
        if(colA~=size(B,2))
            error('MESS:error_arguments',['number of columns of A ' ...
                                'differs with number of columns of B']);
        end
        C = -[(p*(eqn.K_*B(:,one)') + (eqn.K_*B(:,two)'));...
              (eqn.K_*B(:,one)') - p*(eqn.M_*B(:,two)') + eqn.E_*B(:,two)'];
end


end
