function [result, eqn, opts, oper] = init_so_2(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_so_2(eqn, opts, oper, flag1, flag2)
%
% The second order system
%
%   M x"(t) + E x'(t) + K x(t) = B u(t)
%                         y(t) = C x(t)
%       
% is transformed to the first order system
%
%   E_f z'(t) = A_f z(t) + B_f u(t)
%      y(t) = C_f z(t)
%  
% where
%
%        | D  M|
%   E_f= | M  0|
%   
%        |-K  0|
%   A_f= | 0  M|
%   
%         | B |
%   B_f = | 0 |
%   
%   C_f = [C  0]
%   
%         | x(t)  |
%   z(t)= | x'(t) | .
%   
% Matrices M, D, K are assumed to be sparsem quadratic, symmetric
% and positive definit. 
% The function returns true or false if data for A_f and E_f
% respectivey flag1 and flag2  are availabe and correct in the structure eqn. 
%
%   Input:
%
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains function handles for operation with A and E
%   flag1           'A'/'E' to check if A or E is in eqn
%   flag2           'A'/'E' to check if A or E is in eqn
%
%   Output:
%
%   result             1 if data corresponding to flag1 (and flag2)
%                     are available,
%                   0 data are not available  
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains function handles for operation with A and E
%
%   This function does not use other so_2 functions.
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

%start checking
if(nargin<=3)
    error('MESS:control_data','Number of input Arguments are at least 3');

%result = init_so_1(eqn, flag1);    
elseif(nargin==4)
    switch flag1
        case {'A','a'}
            [eqn,result] = checkA(eqn);
        case {'E','e'}
            [eqn,result] = checkE(eqn);
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end
    
%result = init_so_1(eqn,flag1,flag2);
elseif(nargin==5)
    switch flag1
        case {'A','a'}
            [eqn,result] = checkA(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,resultA] = checkA(eqn);
                    result = result && resultA;
                case {'E','e'}
                    [eqn,resultE] = checkE(eqn);
                    result = result && resultE;
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        case {'E','e'}
            [eqn,result] = checkE(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,resultA] = checkA(eqn);
                    result = result &&resultA;
                case {'E','e'}
                    [eqn,resultE] = checkE(eqn);
                    result =  result && resultE;
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end 
end
end

%checkdata for A
function [eqn,result] = checkA(eqn)
result = isfield(eqn,'K_') && isfield(eqn,'E_');
if(result)
    result = isnumeric(eqn.K_) && isnumeric(eqn.E_)...
          && issymmetric(eqn.M_) && issymmetric(eqn.K_);
    if(not(issparse(eqn.K_)) || not(issymmetric(eqn.K_)))
        warning('MESS:control_data','K must be sparse and symmetric');
    end
    if(not(issparse(eqn.E_)) || not(issymmetric(eqn.E_)))
        warning('MESS:control_data','E must be sparse and symmetric');
    end
end
end

%checkdata for E
function [eqn,result] = checkE(eqn)
result = isfield(eqn,'M_')&&isfield(eqn,'K_');
if(result)
    result = isnumeric(eqn.M_) && isnumeric(eqn.K_) ...
          &&issymmetric(eqn.M_) && issymmetric(eqn.K_);
    if(not(issparse(eqn.M_)) || not(issymmetric(eqn.M_)))
        warning('MESS:control_data','M must be sparse and symmetric');
    end
    if(not(issparse(eqn.K_)) || not(issymmetric(eqn.K_)))
        warning('MESS:control_data','K must be sparse and symmetric');
    end
end
end
