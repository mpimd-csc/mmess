function [result, eqn, opts, oper] = ...
    init_state_space_transformed_default(eqn, opts, oper, flag1, flag2)
%% function [result, eqn, opts, oper] = ...
%    init_state_space_transformed_default(eqn, opts, oper, flag1, flag2)
%
% The function returns true or false if data for A_ and E_ resp. flag1 and
% flag2  are availabe and corrects in structure eqn.
%
% Input
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   oper            struct contains function handles for operation
%                   with A and E
%
%   flag1           'A' or 'E' to check if A or E is in eqn
%
%   flag2           'A' or 'E' to check if A or E is in eqn
%
% Output
%   result             1 if data corresponding to flag1 (and flag2) are
%                   available, 0 data are not available
%
%   eqn             structure with data for equations
%
%   opts            structure containing parameter for the algorithm
%
%   oper            struct contains function handles for operation with
%                   A and E
%
% This function does not use other default functions.
%
% This function calls two other functions checkA and checkE implemented
% at the end.
%
% The function checkA(eqn) proofs if a field 'A_' is included in the
% structure eqn and if the field 'A_' is numeric and quadratic.
% This function returns the changed structure eqn and a boolean value
% result (true if 'A_' is in structure eqn and a numeric and quadratic field).
%
% The function checkE(eqn) proofs if a field 'E_' is included in the
% structure eqn and if the field 'E_' is numeric and quadratic.
% If the structure does not include a field E, a new field 'E_' is defined
% as a sparse identity matrix by size of field 'A_'.
% This function returns the changed structure eqn and a boolean value
% result (true if 'E_' is in structure eqn and a numeric and quadratic field).

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%start checking
na = nargin;
if na <= 3
    error( ...
        'MESS:control_data', ...
        'Number of input Arguments are at least 3');

elseif na == 4
    switch flag1
        case {'A', 'a'}
            [eqn, result] = checkA(eqn);
        case {'E', 'e'}
            [eqn, result] = checkE(eqn);
        otherwise
            error( ...
                'MESS:control_data', ...
                'flag1 has to be ''A'' or ''E''');
    end

elseif na == 5
    switch flag1
        case {'A', 'a'}
            [eqn, result] = checkA(eqn);
            switch flag2
                case {'A', 'a'}
                    [eqn, resultA] = checkA(eqn);
                    result = result && resultA;
                case {'E', 'e'}
                    [eqn, resultE]= checkE(eqn);
                    result = result && resultE;
                otherwise
                    error( ...
                        'MESS:control_data', ...
                        'flag2 has to be ''A'' or ''E''');
            end
        case {'E', 'e'}
            [eqn, result] = checkE(eqn);
            switch flag2
                case {'A', 'a'}
                    [eqn, resultA] = checkA(eqn);
                    result = result && resultA;
                case {'E', 'e'}
                    [eqn, resultE]= checkE(eqn);
                    result = result && resultE;
                otherwise
                    error( ...
                        'MESS:control_data', ...
                        'flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error( ...
                'MESS:control_data', ...
                'flag1 has to be ''A'' or ''E''');
    end
end

end

%% Check data for A_.
function [eqn, result] = checkA(eqn)
result = isfield(eqn, 'A_');

if result
    result = isnumeric(eqn.A_);
end

result = result && (size(eqn.A_, 1) == size(eqn.A_, 2));
end

%% Check data for E_.
function [eqn, result] = checkE(eqn)
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end

if not(eqn.haveE)
    result = 1;
    eqn.E_= speye(size(eqn.A_, 1)); % Make sure we have an identity for
                                    % computations in ApE functions.
else
    result = isfield(eqn, 'E_');

    if result
        result = isnumeric(eqn.E_);
    end

    result = result && (size(eqn.E_, 1) == size(eqn.E_, 2));
end

end
