function C = mul_N_default(eqn, opts, opN, B, opB, h) %#ok<INUSL>
% function C = mul_N_default(eqn, opts, opN, B, opB)
%
% This function returns C = N{h}*B, with a given Matrix N{h} and input
% matrix B could be transposed. Matrix N{h} is assumed to be quadratic.
%
%   Inputs:
%
%   eqn     with eqn.N_{h} (cell array) for h = 1,2,...
%
%   opts    structure containing parameters for the algorithm
%
%   opN     character specifying the shape of N{h}
%           opN = 'N' performs N{h} * opB(B)
%           opN = 'T' performs N{h}' * opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifying the shape of B
%           opB = 'N' performs opN(N{h}) * B
%           opB = 'T' performs opN(N{h}) * B'
%
%   h       index h of the current current N{h}
%
%
%   Output:
%
%   C = opN(N) * opB(B)
%
% This function uses another default function size_default(eqn,
% opts) to obtain the number of rows of matrix N in structure eqn.

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check data in eqn structure

rowN = size(eqn.N_{h}, 1);
colN = rowN;

%% perform multiplication
switch opN

    case 'N'
        switch opB

            %implement operation N*B
            case 'N'
                if(colN~=size(B, 1))
                    error('MESS:error_arguments',['number of columns of N ' ...
                        'differs with number of rows of B']);
                end
                C = eqn.N_{h} * B;

                %implement operation N*B'
            case 'T'
                if(colN ~= size(B, 2))
                    error('MESS:error_arguments',['number of columns of N ' ...
                        'differs with number of columns of B']);
                end
                C = eqn.N_{h} * B';
        end

    case 'T'
        switch opB

            %implement operation N'*B
            case 'N'
                if(rowN ~= size(B, 1))
                    error('MESS:error_arguments',['number of rows of N ' ...
                        'differs with number rows of B']);
                end
                C = eqn.N_{h}' * B;

                %implement operatio N'*B'
            case 'T'
                if(rowN ~= size(B, 2))
                    error('MESS:error_arguments',['number of rows of N ' ...
                        'differs with number of columns of B']);
                end
                C = eqn.N_{h}' * B';
        end

end

end

