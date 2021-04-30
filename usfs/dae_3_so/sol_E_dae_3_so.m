function X = sol_E_dae_3_so(eqn, opts, opE, B, opB)%#ok<INUSL>
%% function sol_E_dae_3_so solves opE(E)*X = opB(B) resp. performs X=opE(E)\opB(B)
%
% Input:
%   eqn     structure contains data for E (M_)
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' solves E *X = opB(B)
%           opE = 'T' sovles E'*X = opB(B)
%
%   B       p-x-q matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
% Output
%
%   X       matrix fullfills equation opE(E)*X = opB(B)


%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input Paramters
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
if(not(isfield(eqn, 'M_'))) || not(isnumeric(eqn.M_))
    error('MESS:error_arguments', 'field eqn.M_ is not defined');
end

nv = size(eqn.M_,1);
np = size(eqn.G_,1);

%% solve
if (opB=='N' && (size(B,1)==(2*nv+np))) || (opB=='T' && (size(B,2)==(2*nv+np)))

    switch opE

        case 'N'
            switch opB

                %implement solve E*X=B
                case 'N'
                    X = [B(1:nv,:);
                        [eqn.M_,eqn.alpha*eqn.G_'; ...
                        eqn.alpha*eqn.G_,sparse(np,np)] \ B(nv+1:end,:)
                        ];

                    %implement solve A*X=B'
                case 'T'
                    X = [B(:,1:nv)';
                        [eqn.M_,eqn.alpha*eqn.G_'; ...
                        eqn.alpha*eqn.G_,sparse(np,np)] \ B(:,nv+1:end)'
                        ];

            end

        case 'T'
            switch opB

                %implement solve E'*X=B
                case 'N'
                    X = [B(1:nv,:);
                        [eqn.M_',eqn.alpha*eqn.G_'; ...
                        eqn.alpha*eqn.G_,sparse(np,np)] \ B(nv+1:end,:)
                        ];

                    %implement solve A_'*X=B'
                case 'T'
                    X = [B(:,1:nv)';
                        [eqn.M_',eqn.alpha*eqn.G_'; ...
                        eqn.alpha*eqn.G_,sparse(np,np)] \ B(:,nv+1:end)'
                        ];
            end

    end

elseif (opB=='N' && (size(B,1)==(2*nv))) || (opB=='T' && (size(B,2)==(2*nv)))
    error('MESS:error_usage','sol_E_dae_2_so is only coded for shift parameter computation');
else
    error('MESS:error_arguemnts', 'B has wrong number of cols');
end


end
