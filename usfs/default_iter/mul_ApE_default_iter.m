function C = mul_ApE_default_iter(eqn, opts, opA, p, opE, B, opB)

% function C=mul_ApE_default_iter(eqn, opts,opA,p,opE,B,opB)
%
% This function returns C = (A_+p*E_)*B, where matrices A_ and E_
% given by a structure eqn and input matrix B could be transposed.
%
%   Inputs:
%
%   eqn     structure containing fields 'A_' and 'E_'
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A_
%           opA = 'N' performs (A_ + p*opE(E_))*opB(B)
%           opA = 'T' performs (A_' + p*opE(E_))*opB(B)
%   p       scalar value
%   opE     character specifying the shape of E_
%           opE = 'N' performs (opA(A_) + p*E_)*opB(B)
%           opE = 'T' performs (opA(A_) + p*E_')*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs (opA(A_) + p*opE(E_))*B
%           opB = 'T' performs (opA(A_) + p*opE(E_))*B'
%   Output:
%
%   C = (opA(A_)+ p * opE(E_))*opB(B)
%
% This function uses another default function
% size_default_iter(eqn,opts) to obtain the number of rows of matrix A_
% in structure eqn, that should be equal to the number of rows of
% the matrix E_.
% This function also uses another default function
% mul_A_default_iter(eqn,opA,B,opB) to obtain the result if E=I.
%

% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
%% Check input parameters
if not(ischar(opA)) || not(ischar(opB)) || not(ischar(opE))
    mess_err(opts, 'error_arguments', 'opA, opB or opE is not a char');
end

opA = upper(opA);
opB = upper(opB);
opE = upper(opE);

if not(opA == 'N' || opA == 'T')
    mess_err(opts, 'error_arguments', 'opA is not ''N'' or ''T''');
end

if not(opB == 'N' || opB == 'T')
    mess_err(opts, 'error_arguments', 'opB is not ''N'' or ''T''');
end

if not(opE == 'N' || opE == 'T')
    mess_err(opts, 'error_arguments', 'opE is not ''N'' or ''T''');
end

if (not(isnumeric(p))) || not(length(p) == 1)
    mess_err(opts, 'error_arguments', 'p is not numeric or a scalar value');
end

if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    mess_err(opts, 'error_arguments', 'B has to be a matrix');
end

%% Check data in eqn structure
if eqn.haveE
    if not(isfield(eqn, 'E_')) || not(isfield(eqn, 'A_'))
        mess_err(opts, 'error_arguments', 'field eqn.E_ or eqn.A_ is not defined');
    end
else
    if not(isfield(eqn, 'A_'))
        mess_err(opts, 'error_arguments', 'field eqn.A_ is not defined');
    end
end

rowA = size_default_iter(eqn, opts);
colA = rowA;

if eqn.haveE
    %% Perform multiplication when E_ is not the Identity
    switch opA

        case 'N'
            switch opE

                case 'N'

                    switch opB

                        % Implement operation (A_+p*E_)*B=C
                        case 'N'

                            if not(colA == size(B, 1))
                                mess_err(opts, 'error_arguments', ...
                                         ['number of columns of A_ ' ...
                                          'differs from number of rows of B']);
                            end

                            C = (eqn.A_ + p * eqn.E_) * B;

                            % Implement operation (A_+p*E_)*B'=C
                        case 'T'

                            if not(colA == size(B, 2))
                                mess_err(opts, 'error_arguments', ...
                                         ['number of columns of A_ ' ...
                                          'differs from number of ' ...
                                          'columns of B']);
                            end

                            C = (eqn.A_ + p * eqn.E_) * B';

                    end

                case 'T'

                    switch opB

                        % Implement operation (A_+p*E_')*B=C
                        case 'N'

                            if not(colA == size(B, 1))
                                mess_err(opts, 'error_arguments', ...
                                         ['number of columns ' ...
                                          'of A_ differs from number ' ...
                                          'of rows of B']);
                            end

                            C = (eqn.A_ + p * eqn.E_') * B;

                            % Implement operation (A_+p*E_')*B'=C
                        case 'T'

                            if not(colA == size(B, 2))
                                mess_err(opts, 'error_arguments', ...
                                         ['number of columns ' ...
                                          'of A_ differs from number ' ...
                                          'of columns of B']);
                            end

                            C = (eqn.A_ + p * eqn.E_') * B';

                    end

            end

        case 'T'
            switch opE

                case 'N'

                    switch opB

                        % Implement operation (A_'+p*E_)*B=C
                        case 'N'

                            if not(rowA == size(B, 1))
                                mess_err(opts, 'error_arguments', ...
                                         ['number of rows ' ...
                                          'of A_ differs from number ' ...
                                          'of rows of B']);
                            end

                            C = (eqn.A_' + p * eqn.E_) * B;

                            % Implement operation (A_'+p*E_)*B'=C
                        case 'T'

                            if not(rowA == size(B, 2))
                                mess_err(opts, 'error_arguments', ...
                                         ['number of rows ' ...
                                          'of A_ differs from number ' ...
                                          'of columns of B']);
                            end

                            C = (eqn.A_' + p * eqn.E_) * B';

                    end

                case 'T'

                    switch opB

                        % Implement operation (A_'+p*E_')*B=C
                        case 'N'

                            if not(rowA == size(B, 1))
                                mess_err(opts, 'error_arguments', ...
                                         ['number of rows ' ...
                                          'of A_ differs from number ' ...
                                          'of rows of B']);
                            end

                            C = (eqn.A_' + p * eqn.E_') * B;

                            % Implement operation (A_'+p*E_')*B'=C
                        case 'T'

                            if not(rowA == size(B, 2))
                                mess_err(opts, 'error_arguments', ...
                                         ['number of rows ' ...
                                          'of A_ differs from number ' ...
                                          'of columns of B']);
                            end

                            C = (eqn.A_' + p * eqn.E_') * B';

                    end

            end
    end
elseif not(eqn.haveE)
    %% Perform multiplication when E_ = Identity
    switch opA

        case 'N'

            switch opB

                % Implement operation (A_+p*I)*B=C
                case 'N'

                    if not(colA == size(B, 1))
                        mess_err(opts, 'error_arguments', ...
                                 ['number of columns of ' ...
                                  'A_ differs from number of rows of B']);
                    end

                    C = mul_A_default_iter(eqn, opts, 'N', B, 'N') + p * B;

                    % Implement operation (A_+p*I)*B'=C
                case 'T'

                    if not(colA == size(B, 2))
                        mess_err(opts, 'error_arguments', ...
                                 ['number of columns of ' ...
                                  'A_ differs from number of ' ...
                                  'columns of B']);
                    end

                    C = mul_A_default_iter(eqn, opts, 'N', B, 'T') + p * B';

            end

        case 'T'

            switch opB

                % Implement operation (A_'+p*I)*B=C
                case 'N'

                    if not(rowA == size(B, 1))
                        mess_err(opts, 'error_arguments', ...
                                 ['number of rows of A_ ' ...
                                  'differs from number of rows of B']);
                    end

                    C = mul_A_default_iter(eqn, opts, 'T', B, 'N') + p * B;

                    % Implement operation (A_'+p*I)*B'=C
                case 'T'

                    if not(rowA == size(B, 2))
                        mess_err(opts, 'error_arguments', ...
                                 ['number of rows of A_ ' ...
                                  'differs from number of columns of B']);
                    end

                    C = mul_A_default_iter(eqn, opts, 'T', B, 'T') + p * B';

            end

    end
end
end
