function X = sol_A_dae_2(eqn, opts, opA, B, opB)%#ok<INUSL>
%% function sol_A solves solves opA(A_)*X = opB(B)
% 
% Depending on the vertical dimension of B this solves either with 
%     A = [A1 F;
%           G 0 ]
% or
%  A = P*A1*P'  with
%  P = I - F ( G E1\F ) \ G / E the hidden manifold projector
%  and E1 the corresponding 1,1 block in eqn.E_
%
% Input:
%  eqn       structure with fields A_ and E_
%  opts      struct contains parameters for the algorithm
%  opA       character specifies the form of opA(A_)
%                  opA = 'N' solves A*X=opB(B) 
%                  opA = 'T' solves A^T*X=opB(B)
%
%  B         p-x-q matrix
%
%  opB       character specifies the form of opB(B)
%                  opB = 'N' solves A*X=B 
%                  opB = 'T' solves A*X=B^T
%
% Output:
%  X       matrix fullfills equation opA(A)X = opB(B)
%

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009-2019
%


%% check input Paramters
if (not(ischar(opA)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opA or opB is not a char');
end

opA = upper(opA); opB = upper(opB);
if(not((opA == 'N' || opA == 'T')))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(not(isfield(eqn, 'A_')))
    error('MESS:error_arguments', 'field eqn.A_ is not defined');
end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end

n = size(eqn.A_,1);
st = eqn.st;

[rowB,colB] = size(B);

if(opB == 'N')
    switch rowB
        case n
            dim = n;
        case st
            dim = st;
        otherwise
            error('MESS:error_arguments', 'B has wrong number of rows.');
    end
else
    switch colB
        case n
            dim = n;
        case st
            dim = st;
        otherwise
            error('MESS:error_arguments', 'B has wrong number of columns.');
    end
end

%% solve
if dim ==n
    switch opA
        
        case 'N'
            switch opB
                
                %implement solve A_*X=B
                case 'N'
                    if(n ~= size(B, 1))
                        error('MESS:error_arguments','number of rows of A_ differs with rows of B');
                    end
                    X = eqn.A_ \ B;
                    
                    %implement solve A_*X=B'
                case 'T'
                    if(n ~= size(B, 2))
                        error('MESS:error_arguments','number of rows of A_ differs with cols of B');
                    end
                    X = eqn.A_ \ B';
            end
            
        case 'T'
            switch opB
                
                %implement solve A_'*X=B
                case 'N'
                    if(n ~= size(B, 1))
                        error('MESS:error_arguments','number of cols of A_ differs with rows of B');
                    end
                    X = eqn.A_' \ B;
                    
                    %implement solve A_'*X=B'
                case 'T'
                    if(n ~= size(B, 2))
                        error('MESS:error_arguments','number of cols of A_ differs with cols of B');
                    end
                    X = eqn.A_' \ B';
            end
            
    end
else
     error('MESS:error_arguments','A is singular in these coordinates');
end
