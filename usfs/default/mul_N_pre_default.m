function [eqn, opts, oper] = mul_N_pre_default(eqn, opts, oper)
% Transforms eqn.N_ into a Cell if given as Matrix and counts calls
%
% Input/Output:
%    eqn    struct contains data for equations
%
%    opts   struct contains parameters for the algorithm
%
%    oper   struct contains function handles for operation with N
%
%
% input        eqn.N_           (as matrix or cell)
%
% output       eqn.N_           (as cell)
%              eqn.originalN   (saves matrix version for post_N)
%              eqn.Ncount      (1 for given cell, 0 for given matrix)
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if not(isfield(eqn, 'N_')) || isempty(eqn.N_)
    mess_warn(opts, 'control_data', 'eqn.N_ is missing');
end

% transforms matrix into a cell array
if not(iscell(eqn.N_))

    eqn.originalN = eqn.N_; % saves Matrix for post_N
    rowN = size(eqn.N_, 1);
    colN = size(eqn.N_, 2);

    if rowN == colN
        N = cell (1, 1);
        N{1} = eqn.N_;
    else
        dimN = sqrt(rowN);
        N = cell (1, colN);
        for h = 1:colN
            N{h} = reshape(eqn.N_(:, h), dimN, dimN);
        end
    end

    eqn.Ncount = 1;  % sets flag for post_N (Input was Matrix)
    eqn.N_ = N;

    % no transformation and counts function calls
else
    if not(isfield(eqn, 'Ncount'))
        eqn.Ncount = 1;
    end
    eqn.Ncount = eqn.Ncount + 1;  % sets flag for post_N (Input was Cell)

end

% check if N{h} is quadratic
k = length(eqn.N_);

for h = 1:k
    if not(size(eqn.N_{k}, 1) == size(eqn.N_{k}, 2))
        mess_err(opts, 'error_arguments', ['number of columns of a N{h} ' ...
                                           'differs with number of rows of N{h}']);
    end
end
