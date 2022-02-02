function C = mul_A_dae_1(eqn, opts, opA, B, opB)%#ok<INUSL>
%% function mul_A performs operation C = opA(A_)*opB(B)
%
% Input:
%   eqn     structure contains field A_
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
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%     A = [J1 J2;
%          J3 J4]   J4 regular
%
% Output:
% C = opA(A_)*opB(B)
%
%   uses size_dae_1

%% check input Parameters
if (not(ischar(opA)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opA or opB is not a char');
end

opA = upper(opA); opB = upper(opB);
if(not((opA == 'N' || opA == 'T')))
    error('MESS:error_arguments', 'opA is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to be a matrix');
end

%% check data in eqn structure
if(not(isfield(eqn, 'A_'))) || not(isnumeric(eqn.A_))
    error('MESS:error_arguments', 'field eqn.A_ is not defined');
end

if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end

n = size(eqn.A_,1);
st = eqn.st;
one = 1:st;
two = st + 1 : n;
%% perform multiplication
switch opA

    case 'N'
        switch opB

            %implement operation A_*B
            case 'N'
                if(st > size(B,1))
                    error('MESS:error_arguments', ...
                          'number of cols of A_ differs with rows of B');
                end
                C = eqn.A_(one, one) * B - eqn.A_(one, two) * ...
                    (eqn.A_(two, two) \ (eqn.A_(two, one) * B));

            %implement operation A_*B'
            case 'T'
                if(st > size(B, 2))
                    error('MESS:error_arguments', ...
                          'number of cols of A_ differs with cols of B');
                end
                C = eqn.A_(one, one) * B' - eqn.A_(one, two) * ...
                    (eqn.A_(two, two) \ (eqn.A_(two, one) * B'));
        end

    case 'T'
        switch opB

            %implement operation A_'*B
            case 'N'
                if(st > size(B, 1))
                    error('MESS:error_arguments', ...
                          'number of rows of A_ differs with rows of B');
                end
                C = eqn.A_(one, one)' * B - eqn.A_(two, one)' * ...
                    (eqn.A_(two, two)' \ (eqn.A_(one, two)' * B));

            %implement operatio A_'*B'
            case 'T'
                if(st > size(B, 2))
                    error('MESS:error_arguments', ...
                          'number of rows of A_ differs with cols of B');
                end
                C = eqn.A_(one, one)' * B' - ...
                    eqn.A_(two, one)' * (eqn.A_(two, two)' \ ...
                    (eqn.A_(one, two)' * B'));
        end

end
end
