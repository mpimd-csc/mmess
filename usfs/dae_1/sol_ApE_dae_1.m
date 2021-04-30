function X = sol_ApE_dae_1(eqn, opts, opA, p, opE, B, opB)%#ok<INUSL>

%% function sol_ApE solves
%                  (opA(A_) + p*opE(E_))*X = opB(B)
% resp. performs X=(opA(A_)+p*opE(E_))\opB(B)
%
%
% A_ and E_ are assumed to be quadratic.
% Input:
%
%   eqn     structure contains A_ and E_
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' for A_
%           opA = 'T' for A_'
%
%   p       scalar Value
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' for E_
%           opE = 'T' for E_'
%
%   B       n-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' for B
%           opB = 'T' for B'
%
% Output
%
%   X       matrix fullfills equation (opA(A_)+p*opE(E_))*X = B
%
%   uses size_dae_1

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input Paramters
if (not(ischar(opA)) || not(ischar(opE)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opA, opE or opB is not a char');
end

opA = upper(opA); opE = upper(opE); opB = upper(opB);

if(not((opA == 'N' || opA == 'T')))
    error('MESS:error_arguments', 'opA is not ''N'' or ''T''');
end

if(not((opE == 'N' || opE == 'T')))
    error('MESS:error_arguments', 'opE is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end

if(not(isnumeric(p)))
   error('MESS:error_arguments','p is not numeric');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
if(eqn.haveE)
    if(not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_))...
            || not(isfield(eqn,'A_'))) || not(isnumeric(eqn.A_))
        error('MESS:error_arguments','field eqn.E_ or eqn.A_ is not defined or corrupted');
    end
else
    if(not(isfield(eqn,'A_'))) || not(isnumeric(eqn.A_))
        error('MESS:error_arguments','field eqn.A_ is not defined');
    end
end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end

n = size(eqn.A_,1);
st = eqn.st;
[rowB,colB] = size(B);

if(opB == 'N')
  if (rowB ~= st)
    error('MESS:error_arguments', 'B has not same number of rows as A');
  end
  B = [B;
    zeros(n - st, colB)];
else
  if (colB ~= st)
    error('MESS:error_arguments', 'B has not same number of rows as A');
  end
  B = [B, zeros(rowB, n - st)];
end



%% perfom solve operations for E_ = [ E 0 ]
%                                   [ 0 0 ]
if(eqn.haveE == 1)
    switch opA

        case 'N'
            switch opE

                case 'N'

                    switch opB

                        %implement solve (A_+p*E_)*X=B
                        case 'N'
                            X = (eqn.A_ + p * eqn.E_) \ B;

                        %implement solve (A_+p*E_)*X=B'
                        case 'T'
                            X = (eqn.A_ + p * eqn.E_) \ B';

                    end

                case 'T'

                    switch opB

                        %implement solve (A_+p*E_')*X=B
                        case 'N'
                            X = (eqn.A_ + p * eqn.E_') \ B;

                        %implement solve (A_+p*E_')*X=B'
                        case 'T'
                            X = (eqn.A_ + p * eqn.E_') \ B';

                    end

            end

        case 'T'
            switch opE

                case 'N'

                    switch opB

                        %implement solve (A_'+p*E_)*X=B
                        case 'N'
                            X = (eqn.A_' + p * eqn.E_) \ B;

                        %implement solve (A_'+p*E_)*X=B'
                        case 'T'
                            X = (eqn.A_' + p * eqn.E_) \ B';

                    end

                case 'T'

                    switch opB

                        %implement solve (A_'+p*E_')*X=B
                        case 'N'
                            X = (eqn.A_' + p * eqn.E_') \ B;

                        %implement solve (A_'+p*E_')*X=B'
                        case 'T'
                            X = (eqn.A_' + p * eqn.E_') \ B';

                    end
            end

    end
elseif(eqn.haveE == 0)
    %% perform solve operations for E_ = [ I 0 ]
    %                                    [ 0 0 ]
    switch opA

        case 'N'

            switch opB

                %implement solve (A_+p*E_)*X=B
                case 'N'
                    X = (eqn.A_ + p * eqn.E_) \ B;

                %implement solve (A_+p*E_)*X=B'
                case 'T'
                    X = (eqn.A_ + p * eqn.E_) \ B';

            end

        case 'T'

            switch opB

                %implement solve (A_'+p*E_)*X=B
                case 'N'
                    X = (eqn.A_' + p * eqn.E_) \ B;

                %implement solve (A_'+p*E_)*X=B'
                case 'T'
                    X = (eqn.A_' + p * eqn.E_) \ B';

            end

    end
end
X = X(1 : st, :);
end

