function C = mul_ApE_dae_2(eqn, opts, opA,p,opE, B, opB)%#ok<INUSL>

%% function mul_ApE_default perfoms operation C = (opA(A_)+p*opE(E_))*opB(B)
%
% Input:
%   eqn     structure contains A_
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' performs A_*opB(B)
%           opA = 'T' performs A_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opA(A_)*B
%           opB = 'T' performs opA(A_)*B'
%
% Output:
% C = (opA(A_)+p*opE(E_))*opB(B)
%
%   uses size_dae_2

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%



%% check input Paramters
if (not(ischar(opA)) || not(ischar(opB))|| not(ischar(opE)))
    error('MESS:error_arguments', 'opA, opE or opB is not a char');
end

opA = upper(opA); opB = upper(opB);opE = upper(opE);
if(not((opA == 'N' || opA == 'T')))
    error('MESS:error_arguments', 'opA is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end

if(not((opE == 'N' || opE == 'T')))
    error('MESS:error_arguments', 'opE is not ''N'' or ''T''');
end

if(not(isnumeric(p)))
   error('MESS:error_arguments','p is not numeric');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(not(isfield(eqn, 'A_'))) || not(isnumeric(eqn.A_))
    error('MESS:error_arguments', 'field eqn.A_ is not defined');
end

if(not(isfield(eqn, 'E_'))) || not(isnumeric(eqn.E_))
    error('MESS:error_arguments', 'field eqn.E_ is not defined');
end
n = size(eqn.A_,1);
st = eqn.st;

[rowB,colB] = size(B);

if(opB == 'N')
    if(n > rowB)
        B = [B; zeros(n - st, colB)];
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

%% perfom multiplication
switch opA

  case 'N'

    switch opB

      case 'N'
        %implement operation (A_+p*E_)*B
        C=(eqn.A_+p*eqn.E_)*B;

      case 'T'
        %implement operation (A_+p*E_)*B'
        C=(eqn.A_+p*eqn.E_)*B';
    end

  case 'T'

    switch opB

      case 'N'
        %implement operation (A_+p*E_)'*B
        C=(eqn.A_+p*eqn.E_)'*B;

      case 'T'
        %implement operatio (A_+p*E_)'*B'
        C=(eqn.A_+p*eqn.E_)'*B';
    end

end
end
