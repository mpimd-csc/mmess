function X = sol_A_dae_2_so(eqn, opts, opA, B, opB)%#ok<INUSL>
%  function sol_A solves solves opA(A_)*X = opB(B)
%
% Input:
%  eqn       structure with field A_
%  opts      struct contains parameters for the algorithm
%  opA       character specifies the form of opA(A_)
%                  opA = 'N' solves A_*X=opB(B)
%                  opA = 'T' solves A_^T*X=opB(B)
%  B         p-x-q matrix
%  opB       character specifies the form of opB(B)
%                  opB = 'N' solves A_*X=B
%                  opB = 'T' solves A_*X=B^T
%
% Output:
%  X       matrix fulfills equation opA(A_)X = opB(B)
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%



%% check input Parameters
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
for mat='EKG'
    if(not(isfield(eqn, sprintf('%c_',mat))) || not(eval(sprintf('isnumeric(eqn.%c_))',mat))))
        error('MESS:error_arguments', 'field eqn.%c_ is not defined',mat);
    end
end

nv = size(eqn.M_,1);
np = size(eqn.G_,1);

%% solve

if (opB=='N' && (size(B,1)==(2*nv+np))) || (opB=='T' && (size(B,2)==(2*nv+np)))
switch opA

    case 'N'
        switch opB

            %implement solve A_*X=B
            case 'N'
                x = [eqn.K_, eqn.G_';eqn.G_,sparse(np,np)]\ [B(nv+1:2*nv,:)-eqn.E_*B(1:nv,:);B(2*nv+1:end,:)];
                X = [x(1:nv,:);B(1:nv,:);x(nv+1:end,:)];

            %implement solve A_*X=B'
            case 'T'
                x = [eqn.K_, eqn.G_';eqn.G_,sparse(np,np)]\ [B(:,nv+1:2*nv)'-eqn.E_*B(:,1:nv)';B(:,2*nv+1:end)'];
                X = [x(1:nv,:);B(:,1:nv)';x(nv+1:end,:)];
        end

    case 'T'
        switch opB

            %implement solve A_'*X=B
            case 'N'
                x = [eqn.K_',eqn.G_';eqn.G_,sparse(np,np)]\[B(1:nv,:);B(2*nv+1:end,:)];
                X = [B(nv+1:2*nv,:)-eqn.E_'*x(1:nv,:);x(1:nv,:);x(nv+1:end,:)];

            %implement solve A_'*X=B'
            case 'T'
                x = [eqn.K_',eqn.G_';eqn.G_,sparse(np,np)]\[B(:,1:nv)';B(:,2*nv+1:end)'];
                X = [B(:,nv+1:2*nv)'-eqn.E_'*x(1:nv,:);x(1:nv,:);x(nv+1:end,:)];
        end

end
elseif (opB=='N' && (size(B,1)==(2*nv))) || (opB=='T' && (size(B,2)==(2*nv)))
    error('MESS:error_usage','mul_A_dae_2_so is only coded for shift parameter computation');
else
    error('MESS:error_arguments', 'B has wrong number of cols');
end



end
