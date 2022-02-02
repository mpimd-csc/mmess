function C = mul_ApE_dae_1(eqn, opts, opA, p, opE, B, opB)%#ok<INUSL>
%% function mul_A performs operation C = (opA(A_)+pc*opE(E_))*opB(B)
%
% Input:
%   eqn     structure contains A_
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' performs A_*opB(B)
%           opA = 'T' performs A_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opA(A_)*B
%           opB = 'T' performs opA(A_)*B'
%
% Output:
% C = (opA(A_)+pc*opE(E_))*opB(B)
%
%   uses size_dae_1

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% check input Parameters
if (not(ischar(opA)) || not(ischar(opB)) || not(ischar(opE)))
    error('MESS:error_arguments', 'opA, opB or opE is not a char');
end

opA = upper(opA); opB = upper(opB); opE = upper(opE);

if(not((opA == 'N' || opA == 'T')))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if(not((opE == 'N' || opE == 'T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not(isnumeric(p)))
    error('MESS:error_arguments','p is not numeric');
end

if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(eqn.haveE ==1)
    if(not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_))...
            || not(isfield(eqn,'A_'))) || not(isnumeric(eqn.A_))
        error('MESS:error_arguments', ...
              'field eqn.E_ or eqn.A_ is not defined or corrupted');
    end
else
    if(not(isfield(eqn,'A_'))) || not(isnumeric(eqn.A_))
        error('MESS:error_arguments', ...
              'field eqn.A_ is not defined');
    end
end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
        'Missing or Corrupted st field detected in equation structure.');
end
n = size(eqn.A_,1);
st = eqn.st;
one = 1:st;
two = st + 1 : n;
%% perform multiplication
switch opA

    case 'N'
        switch opE
            case 'N'
                switch opB

                    %implement operation A_*B
                    case 'N'
                        if(st > size(B,1))
                            error('MESS:error_arguments', ...
                                  'number of cols of A_ differs with rows of B');
                        end
                        C = (eqn.A_(one, one) + p * eqn.E_(one, one)) * B ...
                            - eqn.A_(one, two) * (eqn.A_(two, two) \ ...
                            (eqn.A_(two, one) * B));

                    %implement operation A_*B'
                    case 'T'
                        if(st > size(B, 2))
                            error('MESS:error_arguments', ...
                                  'number of cols of A_ differs with cols of B');
                        end
                        C = (eqn.A_(one, one) + p * eqn.E_(one, one)) * B' ...
                            - eqn.A_(one, two) * (eqn.A_(two, two) \ ...
                             (eqn.A_(two, one) * B'));
                end
            case 'T'
                switch opB

                    %implement operation A_*B
                    case 'N'
                        if(st > size(B,1))
                            error('MESS:error_arguments', ...
                                  'number of cols of A_ differs with rows of B');
                        end
                        C = (eqn.A_(one, one) + p * eqn.E_(one, one)') * B ...
                            - eqn.A_(one, two) * (eqn.A_(two, two) \ ...
                              (eqn.A_(two, one) * B));

                    %implement operation A_*B'
                    case 'T'
                        if(st > size(B, 2))
                            error('MESS:error_arguments', ...
                                  'number of cols of A_ differs with cols of B');
                        end
                        C = (eqn.A_(one, one) + p * eqn.E_(one, one)') * B'...
                            - eqn.A_(one, two) * (eqn.A_(two, two) \ ...
                              (eqn.A_(two, one) * B'));
                end
        end

    case 'T'
        switch opE
            case 'N'
                switch opB

                    %implement operation A_'*B
                    case 'N'
                        if(st > size(B, 1))
                            error('MESS:error_arguments', ...
                                  'number of rows of A_ differs with rows of B');
                        end
                        C = (eqn.A_(one, one)' + p * eqn.E_(one, one)) * B ...
                            - eqn.A_(two, one)' * (eqn.A_(two, two)' \ ...
                            (eqn.A_(one, two)' * B));

                    %implement operatio A_'*B'
                    case 'T'
                        if(st > size(B, 2))
                            error('MESS:error_arguments', ...
                                  'number of rows of A_ differs with cols of B');
                        end
                        C = (eqn.A_(one, one)' + p * eqn.E_(one, one)) * B'...
                            - eqn.A_(two, one)' ...
                            * (eqn.A_(two, two)' \ (eqn.A_(one, two)' ...
                            * B'));
                end
            case 'T'
                switch opB

                    %implement operation A_'*B
                    case 'N'
                        if(st > size(B, 1))
                            error('MESS:error_arguments', ...
                                  'number of rows of A_ differs with rows of B');
                        end
                        C = (eqn.A_(one, one)' + p * eqn.E_(one, one)') * B ...
                            - eqn.A_(two, one)' ...
                            * (eqn.A_(two, two)' \ (eqn.A_(one, two)' ...
                            * B));

                    %implement operatio A_'*B'
                    case 'T'
                        if(st > size(B, 2))
                            error('MESS:error_arguments', ...
                                  'number of rows of A_ differs with cols of B');
                        end
                        C = (eqn.A_(one, one)' + p * eqn.E_(one, one)') * B' ...
                            - eqn.A_(two, one)' ...
                            * (eqn.A_(two, two)' \ (eqn.A_(one, two)' ...
                            * B'));
                end
        end

end
end
