function [result, eqn, opts, oper] = init_so_2(eqn, opts, oper, flag1, flag2)
% function [result, eqn, opts, oper] = init_so_2(eqn, opts, oper, flag1, flag2)
%
% Call help mess_usfs_so_2 to see the description of the second order
% system and its transformed first order system
%
%
% The function returns true or false if data for A_f and E_f
% respectively flag1 and flag2  are available and correct in the structure eqn.
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
%   result          1 if data corresponding to flag1 (and flag2)
%                   are available,
%                   0 data are not available
%   eqn             structure with data
%   opts            structure containing parameter for the algorithm
%   oper            struct contains function handles for operation with A and E
%
%   This function does not use other so_2 functions.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
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
