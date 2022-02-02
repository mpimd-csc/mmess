function X = sol_E_dae_1(eqn, opts, opE, B, opB) %#ok<INUSL>
%% function sol_E_dae_1 solves opE(E)*X = opB(B) resp. performs X=opE(E)\opB(B)
%
% Input:
%   eqn     structure contains data for E (M_,K_)
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' solves E *X = opB(B)
%           opE = 'T' solves E'*X = opB(B)
%
%   B       p-x-q matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
% Output
%
%   X       matrix fulfills equation opE(E)*X = opB(B)
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
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end
st = eqn.st;

%% solve
switch opE

    case 'N'
        switch opB

            %implement solve A_*X=B
            case 'N'
                if(st ~= size(B, 1))
                    error('MESS:error_arguments','number of rows of A_ differs with rows of B');
                end
                X = eqn.E_(1 : st, 1 : st) \ B;

            %implement solve A_*X=B'
            case 'T'
                if(st ~= size(B, 2))
                    error('MESS:error_arguments','number of rows of A_ differs with cols of B');
                end
                X = eqn.E_(1 : st, 1 : st) \ B';
        end

    case 'T'
        switch opB

            %implement solve A_'*X=B
            case 'N'
                if(st ~= size(B, 1))
                    error('MESS:error_arguments','number of cols of A_ differs with rows of B');
                end
                X = eqn.E_(1 : st, 1 : st)' \ B;

            %implement solve A_'*X=B'
            case 'T'
                if(st ~= size(B, 2))
                    error('MESS:error_arguments','number of cols of A_ differs with cols of B');
                end
                X = eqn.E_(1 : st, 1 : st)' \ B';
        end

end

end
