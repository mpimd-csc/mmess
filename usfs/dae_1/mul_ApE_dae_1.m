function C = mul_ApE_dae_1(eqn, opts, opA, p, opE, B, opB)
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
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check input Parameters
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

if not(isnumeric(p))
    mess_err(opts, 'error_arguments', 'p is not numeric');
end

if not(isfield(eqn, 'haveE'))
    eqn.haveE = false;
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    mess_err(opts, 'error_arguments', 'B has to ba a matrix');
end

%% check data in eqn structure
if eqn.haveE
    if (not(isfield(eqn, 'E_')) || not(isnumeric(eqn.E_)) || ...
        not(isfield(eqn, 'A_'))) || not(isnumeric(eqn.A_))
        mess_err(opts, 'error_arguments', ...
                 'field eqn.E_ or eqn.A_ is not defined or corrupted');
    end
else
    if (not(isfield(eqn, 'A_'))) || not(isnumeric(eqn.A_))
        mess_err(opts, 'error_arguments', ...
                 'field eqn.A_ is not defined');
    end
end
if not(isfield(eqn, 'manifold_dim')) || not(isnumeric(eqn.manifold_dim))
    mess_err(opts, 'equation_data', ...
             ['Missing or corrupted manifold_dim field detected in ' ...
              'equation structure.']);
end
n = size(eqn.A_, 1);

one = 1:eqn.manifold_dim;
two = eqn.manifold_dim + 1:n;
%% perform multiplication
switch opA

    case 'N'
        switch opE
            case 'N'
                switch opB

                    % implement operation A_*B
                    case 'N'
                        if eqn.manifold_dim > size(B, 1)
                            mess_err(opts, 'error_arguments', ...
                                     ['number of cols of A_ differs ' ...
                                      'from rows of B']);
                        end
                        C = (eqn.A_(one, one) + p * eqn.E_(one, one)) * B - ...
                             eqn.A_(one, two) * (eqn.A_(two, two) \ ...
                                                 (eqn.A_(two, one) * B));

                        % implement operation A_*B'
                    case 'T'
                        if eqn.manifold_dim > size(B, 2)
                            mess_err(opts, 'error_arguments', ...
                                     ['number of cols of A_ differs ' ...
                                      'from cols of B']);
                        end
                        C = (eqn.A_(one, one) + p * eqn.E_(one, one)) * B' - ...
                            eqn.A_(one, two) * (eqn.A_(two, two) \ ...
                                                (eqn.A_(two, one) * B'));
                end
            case 'T'
                switch opB

                    % implement operation A_*B
                    case 'N'
                        if eqn.manifold_dim > size(B, 1)
                            mess_err(opts, 'error_arguments', ...
                                     ['number of cols of A_ differs ' ...
                                      'from rows of B']);
                        end
                        C = (eqn.A_(one, one) + p * eqn.E_(one, one)') * B - ...
                            eqn.A_(one, two) * (eqn.A_(two, two) \ ...
                                                (eqn.A_(two, one) * B));

                        % implement operation A_*B'
                    case 'T'
                        if eqn.manifold_dim > size(B, 2)
                            mess_err(opts, 'error_arguments', ...
                                     ['number of cols of A_ differs ' ...
                                      'from cols of B']);
                        end
                        C = (eqn.A_(one, one) + p * eqn.E_(one, one)') * ...
                            B' - eqn.A_(one, two) * (eqn.A_(two, two) \ ...
                                                     (eqn.A_(two, one) * B'));
                end
        end

    case 'T'
        switch opE
            case 'N'
                switch opB

                    % implement operation A_'*B
                    case 'N'
                        if eqn.manifold_dim > size(B, 1)
                            mess_err(opts, 'error_arguments', ...
                                     ['number of rows of A_ differs ' ...
                                      'with rows of B']);
                        end
                        C = (eqn.A_(one, one)' + p * eqn.E_(one, one)) * B - ...
                            eqn.A_(two, one)' * (eqn.A_(two, two)' \ ...
                                                 (eqn.A_(one, two)' * B));

                        % implement operatio A_'*B'
                    case 'T'
                        if eqn.manifold_dim > size(B, 2)
                            mess_err(opts, 'error_arguments', ...
                                     ['number of rows of A_ differs ' ...
                                      'with cols of B']);
                        end
                        C = (eqn.A_(one, one)' + p * eqn.E_(one, one)) * ...
                            B' - eqn.A_(two, one)' * ...
                            (eqn.A_(two, two)' \ (eqn.A_(one, two)' * B'));
                end
            case 'T'
                switch opB

                    % implement operation A_'*B
                    case 'N'
                        if eqn.manifold_dim > size(B, 1)
                            mess_err(opts, 'error_arguments', ...
                                     ['number of rows of A_ differs ' ...
                                      'with rows of B']);
                        end
                        C = (eqn.A_(one, one)' + p * eqn.E_(one, one)') * ...
                            B - eqn.A_(two, one)' * ...
                            (eqn.A_(two, two)' \ (eqn.A_(one, two)' * B));

                        % implement operatio A_'*B'
                    case 'T'
                        if eqn.manifold_dim > size(B, 2)
                            mess_err(opts, 'error_arguments', ...
                                     ['number of rows of A_ differs ' ...
                                      'with cols of B']);
                        end
                        C = (eqn.A_(one, one)' + p * eqn.E_(one, one)') * ...
                            B' - eqn.A_(two, one)' * ...
                            (eqn.A_(two, two)' \ (eqn.A_(one, two)' * B'));
                end
        end

end
end
