function C = mul_E_dae_1(eqn, opts, opE, B, opB)%#ok<INUSL>

%% function mul_A performs operation C = opE(E_)*opB(B)
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
%   uses no other dae_1 function

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input Parameters
if not(ischar(opE)) || not(ischar(opB))
    error('MESS:error_arguments', 'opE or opB is not a char');
end

opE = upper(opE); opB = upper(opB);
if not(opE == 'N' || opE == 'T')
    error('MESS:error_arguments', 'opE is not ''N'' or ''T''');
end

if not(opB == 'N' || opB == 'T')
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end

if not(isnumeric(B)) || not(ismatrix(B))
    error('MESS:error_arguments', 'B has to ba a matrix');
end
%% check data in eqn structure
if not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_))
    error('MESS:error_arguments', ...
        'Missing or Corrupted E_ field detected in equation structure.');
end
if not(isfield(eqn, 'st')) || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end
st = eqn.st;
one = 1:st;
%% perform multiplication
switch opE

    case 'N'
        switch opB

            %implement operation E_*B
            case 'N'
                if(st ~= size(B, 1))
                    error('MESS:error_arguments','number of cols of E_ differs with rows of B');
                end
                C = eqn.E_(one, one) * B;

            %implement operation E_*B'
            case 'T'
                if(st ~= size(B, 2))
                    error('MESS:error_arguments','number of cols of E_ differs with cols of B');
                end
                C = eqn.E_(one, one) * B';
        end

    case 'T'
        switch opB

            %implement operation E_'*B
            case 'N'
                if(st ~= size(B, 1))
                    error('MESS:error_arguments','number of rows of E_ differs with rows of B');
                end
                C = eqn.E_(one, one)' * B;

            %implement operatio E_'*B'
            case 'T'
                if(st ~= size(B, 2))
                    error('MESS:error_arguments','number of rows of E_ differs with cols of B');
                end
                C = eqn.E_(one, one)' * B';
        end

end

end
