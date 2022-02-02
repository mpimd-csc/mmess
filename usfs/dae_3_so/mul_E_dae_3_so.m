function C = mul_E_dae_3_so(eqn, opts, opE, B, opB)%#ok<INUSL>
%% function mul_A performs operation C = opE(E_)*opB(B)
%
% Input:
%   eqn     structure contains field E_
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' performs E_*opB(B)
%           opE = 'T' performs E_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opE(E_)*B
%           opB = 'T' performs opE(E_)*B'
%
% Output:
% C = opE(E_)*opB(B)
%

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
for mat='MG'
    if(not(isfield(eqn, sprintf('%c_',mat))) || ...
       not(eval(sprintf('isnumeric(eqn.%c_)',mat))))
        error('MESS:error_arguments', 'field eqn.%c_ is not defined',mat);
    end
end

%% perform multiplication
nv = size(eqn.M_,1);
np = size(eqn.G_,1);

if (opB=='N' && (size(B,1)==(2*nv+np))) || (opB=='T' && (size(B,2)==(2*nv+np)))
    switch opE

        case 'N'
            switch opB

                %implement operation E_*B
                case 'N'
                    C = [B(1:nv,:);
                        eqn.M_*B(nv+1:2*nv,:)+eqn.alpha*eqn.G_'*B(2*nv+1:end,:);
                        eqn.alpha*eqn.G_*B(nv+1:2*nv,:)];

                    %implement operation E_*B'
                case 'T'
                    C = [B(:,1:nv)';
                        eqn.M_*B(:,nv+1:2*nv)'+eqn.alpha*eqn.G_'*B(:,2*nv+1:end)';
                        eqn.alpha*eqn.G_*B(:,nv+1:2*nv)'];

            end

        case 'T'

            switch opB
                %implement operation E_'*B
                case 'N'
                    C = [B(1:nv,:);
                        eqn.M_'*B(nv+1:2*nv,:)+eqn.alpha*eqn.G_'*B(2*nv+1:end,:);
                        eqn.alpha*eqn.G_*B(nv+1:2*nv,:)];

                    %implement operation E_'*B'
                case 'T'
                    C = [B(:,1:nv)';
                        eqn.M_'*B(:,nv+1:2*nv)'+eqn.alpha*eqn.G_'*B(:,2*nv+1:end)';
                        eqn.alpha*eqn.G_*B(:,nv+1:2*nv)'];

            end

    end

elseif (opB=='N' && (size(B,1)==(2*nv))) || (opB=='T' && (size(B,2)==(2*nv)))
    switch opE

        case 'N'
            switch opB

                %implement operation E_*B
                case 'N'
                    C = [B(1:nv,:);
                        eqn.M_*B(nv+1:2*nv,:)];

                    %implement operation E_*B'
                case 'T'
                    C = [B(:,1:nv)';
                        eqn.M_*B(:,nv+1:2*nv)'];

            end

        case 'T'

            switch opB
                %implement operation E_'*B
                case 'N'
                    C = [B(1:nv,:);
                        eqn.M_'*B(nv+1:2*nv,:)];

                    %implement operation E_'*B'
                case 'T'
                    C = [B(:,1:nv)';
                        eqn.M_'*B(:,nv+1:2*nv)'];
            end

    end
else
    error('MESS:error_arguments', 'B has wrong number of cols');
end

end
