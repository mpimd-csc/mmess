function C = mul_E_dae_2(eqn, opts, opE, B, opB)%#ok<INUSL>

%% function mul_E performs operation C = opE(E_)*opB(B)
%
% Input:
%   eqn     structure contains field E_
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' performs E_*opB(B)
%           opE = 'T' performs E_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opE(E_)*B
%           opB = 'T' performs opE(E_)*B'
%
% Output:
% C = opE(E_)*opB(B)
%
%   uses no other dae_2 function

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input Parameters
if (not(ischar(opE)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opE or opB is not a char');
end

opE = upper(opE); opB = upper(opB);
if(not((opE == 'N' || opE == 'T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end
%% check data in eqn structure
if(not(isfield(eqn, 'E_'))) || not(isnumeric(eqn.E_))
  error('MESS:error_arguments', 'field eqn.E_ is not defined');
end

if(not(isfield(eqn, 'S_'))) || not(isnumeric(eqn.S_))
    error('MESS:error_arguments', ['field eqn.S_ is not defined. Did ' ...
                        'you forget to run mul_E_pre?']);
end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end

st = eqn.st;
switch opB
  case 'N'
    rowB=size(B,1);
  case 'T'
    rowB=size(B,2);
end
if rowB~=st && rowB~=size(eqn.E_,1)
  error('MESS:error_arguments', 'size of B does not match data in E');
end

%% perform multiplication
switch opE

    case 'N'
        switch opB

            %implement operation E_*B
            case 'N'
              C = eqn.S_(1 : rowB, 1 : rowB) * B;

            %implement operation E_*B'
            case 'T'
              C = eqn.S_(1 : rowB, 1 : rowB) * B';
        end

    case 'T'
        switch opB

            %implement operation E_'*B
            case 'N'
                C = eqn.S_(1 : rowB, 1 : rowB)' * B;

            %implement operation E_'*B'
            case 'T'
                C = eqn.S_(1 : rowB, 1 : rowB)' * B';
        end

end
% This portion would make multiplication with E more correct. Still,
% currently explicit projection is not needed anywhere in our codes and it
% easily double the runtime.
% if rowB==st
%     C = mul_Pi(eqn,'N',C,'N');
% end
end
