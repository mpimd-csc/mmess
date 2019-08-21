function C = mul_E_state_space_transformed_default(eqn, opts, opE, B, opB)
%% function C = mul_E_state_space_transformed_default(eqn,opts,opE,B,opB)
%
% This function returns C = E_*B, where matrix E_ given by structure eqn
% and input matrix B could be transposed.
% Matrix E_ is assumed to be quadratic.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   opE             character specifying the shape of E_
%                       opE = 'N' performs E_*opB(B)
%                       opE = 'T' performs E_'*opB(B)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                       opB = 'N' performs opE(E_)*B
%                       opB = 'T' performs opE(E_)*B'
%
% Output
%   C = opE(E_)*opB(B)
%
% This function uses another default function size_default(eqn, opts) to
% obtain the number of rows of matrix E_ in structure eqn.

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

%% Check input parameters.
assert(ischar(opE) && ischar(opB), ...
    'MESS:error_arguments', ...
    'opE or opB is not a char');

opE = upper(opE);
opB = upper(opB);

assert((opE == 'N') || (opE == 'T'), ...
    'MESS:error_arguments', ...
    'opE is not ''N'' or ''T''');

assert((opB == 'N') || (opB == 'T'), ...
    'MESS:error_arguments', ...
    'opB is not ''N'' or ''T''');

assert(isnumeric(B) && ismatrix(B), ...
    'MESS:error_arguments', ...
    'B has to ba a matrix');

%% Check data in eqn structure.
assert(isfield(eqn,'E_'), ...
    'MESS:error_arguments', ...
    'field eqn.E_ is not defined');

rowE = size_default(eqn, opts);
colE = rowE;

%% Perform multiplication.
switch opE
    case 'N'
        switch opB
            case 'N' % Implement operation E_*B.
                assert(colE == size(B, 1), ...
                    'MESS:error_arguments', ...
                    ['number of columns of E_ differs with ' ...
                    'number of rows of B']);
                C = eqn.E_ * B;
                
            case 'T' % Implement operation E_*B'.
                assert(colE == size(B, 2), ...
                    'MESS:error_arguments', ...
                    ['number of columns of E_ differs with ' ...
                    'number of columns of B']);
                C = eqn.E_ * B';
        end
        
    case 'T'
        switch opB
            case 'N' % Implement operation E_'*B.
                assert(rowE == size(B, 1), ...
                    'MESS:error_arguments', ...
                    ['number of rows of E_ differs with ' ...
                    'number rows of B']);
                C = eqn.E_' * B;
                
            case 'T' % Implement operation E_'*B'.
                assert(rowE == size(B, 2), ...
                    'MESS:error_arguments', ...
                    ['number of rows of E_ differs with ' ...
                    'number of columns of B']);
                C = eqn.E_' * B';
        end
        
end

end
