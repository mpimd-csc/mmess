function n = size_default(eqn, opts, oper)%#ok<INUSD>
% function n = size_default(eqn, opts, oper)
%
% This function returns the number of rows of matrix A_ in structure eqn.
%
%   Input:
%
%   eqn     struct contains data for equations
%
%   opts    struct contains parameters for the algorithm
%
%   oper    struct contains function handles for operation 
%           with A and E
%
%   Output:
%
%   n       number of rows of matrix A_ in structure eqn
%
% This function does not use other default functions.

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
if(not(isfield(eqn,'A_')))
    error('MESS:error_arguments','field eqn.A_ is not defined');
end

n = size(eqn.A_,1);
end

