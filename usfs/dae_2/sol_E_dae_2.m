function X = sol_E_dae_2(eqn, opts, opE, B, opB) %#ok<INUSL>
%% function sol_E solves opE(S_)*X = opB(B) resp. performs X=opE(S_)\opB(B)
% sol_E_pre should be called before to construct
% S_ = [ E1 -J';
%      [ J   0 ]
% from
% A = [ A1 -J';
%       J   0]
% E = [ E1 0;
%       0  0]
%
% Input:
%   eqn     structure contains data for S_
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
%   uses no other dae_2 function

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
if(not(isfield(eqn, 'S_')))
    error('MESS:error_arguments', ['field eqn.S_ is not defined. Did ' ...
                        'you forget to run sol_E_pre?']);
end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end

st = eqn.st;

n = size(eqn.S_,1);

[rowB,colB] = size(B);

if(opB == 'N')
    if(rowB == st)
        B = [B; zeros(n - rowB, colB)];
    elseif rowB ~= n
        error('MESS:error_arguments', 'size of B does not match data in E');
    end
else
    if(colB == st)
        B = [B, zeros(rowB, n - colB)];
    elseif colB ~= n
        error('MESS:error_arguments', 'size of B does not match data in E');
    end
end


%% solve
switch opE

    case 'N'
        switch opB

            %implement solve S_*X=B
            case 'N'

                X = eqn.S_ \ B;

            %implement solve S_*X=B'
            case 'T'

                X = eqn.S_ \ B';
        end

    case 'T'
        switch opB

            %implement solve S_'*X=B
            case 'N'

                X = eqn.S_' \ B;

            %implement solve S_'*X=B'
            case 'T'

                X = eqn.S_' \ B';
        end

end

end
