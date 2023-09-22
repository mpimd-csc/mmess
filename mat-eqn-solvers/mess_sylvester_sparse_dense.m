function [X, eigH, eqn, opts, oper] = mess_sylvester_sparse_dense(varargin)
%% MESS_SYLVESTER_SPARSE_DENSE Solves the Sylvester equation
%
%                A * X * F + E * X * H = -M                        (1)
%
% with F, H small and dense and A, E large and sparse.
%
% Calling sequence:
%
% [X ,eigH, eqn, ops, oper] = ...
%      mess_sylvester_sparse_dense(A, TransA, H, TransH, M, E, F)
%
% or
%
% [X ,eigH, eqn, ops, oper] = ...
%      mess_sylvester_sparse_dense(eqn, opts,oper, H, TransH, M, F)
%
% If TransA == 'T' the matrices A and E are treated transposed in (1).
% Similarly, TransH == 'T' indicates the same for F, H.
%
% The inputs E and F are optional and default to the identities of
% appropriate size when omitted.
%
% If TransA, TransH are not 'T', or not set, they are set to 'N', i.e.,
% (1) is solved just as given above. In the "eqn, opts, oper" case TransA
% is decided from eqn.type and existence of E is determined from eqn.haveE.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check inputs
opts = struct;
if isstruct(varargin{1}) && ...
       isstruct(varargin{2}) && ...
       isstruct(varargin{3})

    % The case: mess_sylvester_sparse_dense(eqn, opts, oper, ...
    %                                       H, TransH, M, F)
    if nargin < 6
        mess_err(opts, 'inputs', ...
                 ['At least 6 arguments required, ', ...
                  'when first input is a matrix']);

    end

    eqn  = varargin{1};
    opts = varargin{2};
    oper = varargin{3};

    if ismatrix(varargin{4})
        H = full(varargin{4});
    else
        mess_err(opts, 'inputs', ...
                 ['Fourth argument must be a dense matrix, ', ...
                  'when first is a matrix']);
    end

    if isa(varargin{5}, 'char')
        TransH = varargin{5};
    else
        mess_err(opts, 'inputs', ...
                 'Fifth argument must be a char, when first is a matrix');
    end
    if not(TransH == 'T')
        TransH = 'N';
    end

    if ismatrix(varargin{6})
        M = varargin{6};
    else
        mess_err(opts, 'inputs', ...
                 ['Sixth argument must be a matrix, ', ...
                  'when the first is a matrix']);
    end

    if nargin < 7
        haveF = false;
        F = [];
    elseif ismatrix(varargin{7})
        haveF = true;
        F = full(varargin{7});
    else
        mess_err(opts, 'inputs', ...
                 ['7th argument must be a matrix, ', ...
                  'when the first is a matrix']);
    end

elseif isnumeric(varargin{1}) && ismatrix(varargin{1})
    % The case: mess_sylvester_sparse_dense(A, TransA, H, TransH, M, E, F)
    if nargin < 5
        mess_err(opts, 'inputs', ...
                 ['At least 5 arguments required, ', ...
                  'when first input is a matrix']);

    end

    if isa(varargin{2}, 'char')
        TransA = varargin{2};
    else
        mess_err(opts, 'inputs', ...
                 'Second argument must be a char, when first is a matrix');
    end
    if not(TransA == 'T')
        TransA = 'N';
    end

    if ismatrix(varargin{3})
        H = full(varargin{3});
    else
        mess_err(opts, 'inputs', ...
                 ['Third argument must be a dense matrix, ', ...
                  'when first is a matrix']);
    end

    if isa(varargin{4}, 'char')
        TransH = varargin{4};
    else
        mess_err(opts, 'inputs', ....
                 'Second argument must be a char, when first is a matrix');
    end
    if not(TransH == 'T')
        TransH = 'N';
    end

    if ismatrix(varargin{5})
        M = varargin{5};
    else
        mess_err(opts, 'inputs', ...
                 ['Fifth argument must be a matrix, ', ...
                  'when the first is a matrix']);
    end

    if nargin < 6
        eqn.haveE = false;
    else
        eqn.E_ = varargin{6};
        % If the user passed the identity for E we still want haveE to be
        % false
        eqn.haveE = not(full(sum(sum(eqn.E_ - speye(size(eqn.E_, 1)) > ...
                                     eps))) == 0);
        if not(eqn.haveE)
            eqn = rmfield(eqn, 'E_');
        end
    end

    if nargin < 7
        haveF = false;
        F = [];
    elseif ismatrix(varargin{7})
        haveF = true;
        F = full(varargin{7});
    else
        mess_err(opts, 'inputs', ...
                 ['7th argument must be a matrix, ', ...
                  'when the first is a matrix']);
    end

    % set eqn
    eqn.A_ = varargin{1};
    eqn.type = TransA;

    % set oper and opts
    [oper, opts] = operatormanager(opts, 'default');
    opts.norm = 2;

end

%% Transformation to triangular form
if not(haveF)
    [Q, S] = schur(H, 'complex');
    M = M * Q;
else
    [T, S, Q, Z] = qz(F, H);
    switch TransH
        case 'N'
            M = M * Z;
            Q = Q';
        case 'T'
            M = M * Q';
            Q = Z;
    end
end

[~, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');
[eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);

nH = size(H, 1);
nA = oper.size(eqn, opts, oper);
X  = zeros(nA, nH);
AX = zeros(nA, nH);
EX = zeros(nA, nH);

if TransH == 'N'
    if not(haveF) % case: AX + EXH + M = 0
        for j = 1:nH
            prev = 1:(j - 1);
            rhs = X(:, prev) * S(prev, j);
            rhs = -M(:, j) - oper.mul_E(eqn, opts, eqn.type, rhs, 'N');
            X(:, j) = oper.sol_ApE(eqn, opts, eqn.type, ...
                                   S(j, j), eqn.type, rhs, 'N');
        end

    else % case: AXF + EXH + M = 0
        for j = 1:nH
            prev = 1:(j - 1);
            rhs = -M(:, j) - ...
                   AX(:, prev) * T(prev, j) - ...
                   EX(:, prev) * S(prev, j);
            X(:, j) = oper.sol_ApE(eqn, opts, ...
                                   eqn.type, ...
                                   S(j, j) / T(j, j), ...
                                   eqn.type, ...
                                   T(j, j) \ rhs, ...
                                   'N');

            if j < nH
                AX(:, j) = oper.mul_A(eqn, opts, eqn.type, X(:, j), 'N');
                EX(:, j) = oper.mul_E(eqn, opts, eqn.type, X(:, j), 'N');
            end
        end
    end
else % TransH == 'T'
    if not(haveF)   % case: A'X + E'XH + M = 0
        for j = nH:-1:1
            prev = (j + 1):nH;
            rhs = X(:, prev) * S(j, prev)';
            rhs = -M(:, j) - oper.mul_E(eqn, opts, eqn.type, rhs, 'N');
            X(:, j) = oper.sol_ApE(eqn, opts, ...
                                   eqn.type, ...
                                   conj(S(j, j)), ...
                                   eqn.type, ...
                                   rhs, ...
                                   'N');
        end

    else    % case: A'XF + E'XH + M = 0
        for j = nH:-1:1
            prev = (j + 1):nH;
            rhs = -M(:, j) - ...
                   AX(:, prev) * T(j, prev)' - ...
                   EX(:, prev) * S(j, prev)';
            X(:, j) = oper.sol_ApE(eqn, opts, ...
                                   eqn.type, ...
                                   conj(S(j, j)) / conj(T(j, j)), ...
                                   eqn.type, ...
                                   conj(T(j, j)) \ rhs, ...
                                   'N');
            if j > 1
                AX(:, j) = oper.mul_A(eqn, opts, eqn.type, X(:, j), 'N');
                EX(:, j) = oper.mul_E(eqn, opts, eqn.type, X(:, j), 'N');
            end
        end
    end
end

X = X * Q';

if isreal(H) && isreal(F)
    X = real(X);
end
eigH = diag(S);

[eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);

end
