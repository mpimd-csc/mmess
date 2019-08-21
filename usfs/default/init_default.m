function [result, eqn, opts, oper] = init_default(eqn, opts, oper, flag1, flag2)

% function [result, eqn, opts, oper] = init_default(eqn, opts, oper, flag1, flag2)
%
% The function returns true or false if data for A_ and E_
% resp. flag1 and flag2  are availabe and corrects in structure
% eqn. 
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
%   are available , 0 data are not available
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains function handles for operation with A and E
%
%   This function does not use other default functions.
%
%   This function calls two other functions checkA and checkE
%   implemented at the end. 
%
%   The function checkA(eqn) proofs if a field 'A_' is included in
%   the structure eqn and if the field 'A_' is numeric and
%   quadratic. 
%   This function returns the changed structure eqn and a boolean
%   value result (1- 'A_' is in structure eqn and a numeric and
%   quadratic field) 
%
%   The function checkE(eqn) proofs if a field 'E_' is included in
%   the structure eqn and if the field 'E_' is numeric and
%   quadratic. 
%   If the structure does not include a field E, a new field 'E_'
%   is defined as a sparse identity matrix by size of field 'A_'. 
%   This function returns the changed structure eqn and a boolean
%   value result (1- 'E_' is in structure eqn and a numeric and
%   quadratic field) 

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
na = nargin;
if not(isfield(eqn, 'LTV')),  eqn.LTV=0; end
if(na<3)
    error('MESS:control_data','Number of input Arguments are at least 4');
    
    %result = init_default(eqn, flag1);    
elseif(nargin==4)
    switch flag1
      case {'A','a'}
            if eqn.LTV
                [eqn, result] = checkA_time(eqn,opts);
            else
                [eqn, result] = checkA(eqn);
            end
      case {'E','e'}
            if eqn.LTV
                [eqn, result] = checkE_time(eqn,opts);
            else
                [eqn, result] = checkE(eqn);
            end
      otherwise
        error('MESS:control_data','flag1 has to be ''A_'' or ''E_''');
    end
    
    %result = init_default(eqn,flag1,flag2);
elseif(nargin==5)
    switch flag1
      case {'A','a'}
            if eqn.LTV
                [eqn, result] = checkA_time(eqn,opts);
            else
                [eqn, result] = checkA(eqn);
            end
        switch flag2
          case {'A','a'}
                    if eqn.LTV
                        [eqn, resultA] = checkA_time(eqn,opts);
                    else
                        [eqn, resultA] = checkA(eqn);
                    end
                    result = result && resultA;
          case {'E','e'}
                    if eqn.LTV
                        [eqn, resultE] = checkE_time(eqn,opts);
                    else
                        [eqn, resultE]= checkE(eqn);
                    end
                    result = result && resultE;
          otherwise
            error('MESS:control_data','flag2 has to be ''A'' or ''E''');
        end
      case {'E','e'}
            if eqn.LTV
                [eqn, result] = checkE_time(eqn,opts);
            else
                [eqn, result] = checkE(eqn);
            end
        switch flag2
          case {'A','a'}
                    if eqn.LTV
                        [eqn, resultA] = checkA_time(eqn,opts);
                    else
                        [eqn, resultA] = checkA(eqn);
                    end
                    result = result && resultA;
          case {'E','e'}
                    if eqn.LTV
                        [eqn,resultE] = checkE_time(eqn,opts);
                    else
                        [eqn, resultE]= checkE(eqn);
                    end
                    result = result && resultE;
          otherwise
            error('MESS:control_data','flag2 has to be ''A'' or ''E''');
        end
      otherwise
        error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end
end
end

%checkdata for A_
function [eqn, result] = checkA(eqn)
result = isfield(eqn,'A_');
if(result)
    result = isnumeric(eqn.A_);
end
result=result&&(size(eqn.A_,1)==size(eqn.A_,2));
end

%checkdata for E_
function [eqn, result] = checkE(eqn)
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end

if not(eqn.haveE)
    result = 1;
    eqn.E_= speye(size(eqn.A_,1)); %make sure we have an identity for
                                   %computations in ApE functions
else
    result = isfield(eqn,'E_');
    if(result)
        result = isnumeric(eqn.E_);
    end
    result=result&&(size(eqn.E_,1)==size(eqn.E_,2));
end
end

%checkdata for A_time
function [eqn, result] = checkA_time(eqn,opts)
result = isa(eqn.A_time,'function_handle');
A = eqn.A_time(opts.t0);
if(result)
    result = isnumeric(A);
end
result=result&&(size(A,1)==size(A,2));
end

%checkdata for E_time
function [eqn, result] = checkE_time(eqn,opts)
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
if not(eqn.haveE)
    result = 1;
    eqn.E_= speye(size(eqn.A_time(opts.t0),1)); %make sure we have an identity for
    %computations in ApE functions
else
    result = isa(eqn.E_time,'function_handle') ...
        && isa(eqn.dt_E_time,'function_handle');
    E = eqn.E_time(opts.t0);
    if(result)
        result = isnumeric(E);
    end
    result=result&&(size(E,1)==size(E,2));
end
end
