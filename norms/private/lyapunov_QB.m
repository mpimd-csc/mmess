function y = lyapunov_QB(Z, x, eqn, oper, opts, ~)
% Computes matrix vector product with the Lyapunov operator.
%
% Input:
%  Z         Low-rank solution factor of the Riccati equation
%
%  x         vector for matrix vector product
%
%  eqn       structure with data for A, E, B and C
%
%  oper      structure contains function handles for operations with
%                  A, E and N
%
%  opts      full options structure (passed on to function handles in oper)
%
%  D         considered empty
%
%  eqn.type   'N' or 'T' for the type of Lyapunov equation
%
% Output:
%  y         result of matrix vector product

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Setting eqn.type
switch eqn.type
    case 'N'
        adjoint = 'T';
    case 'T'
        adjoint = 'N';
    otherwise
        mess_err(opts, 'control_data', 'eqn.type has to be ''N'' or ''T''');
end

%% Z*(Z'*E)
if eqn.haveE
    z = Z * (Z' * (oper.mul_E(eqn, opts, adjoint, x, 'N')));
else
    z = Z * (Z' * x);
end
y1 = oper.mul_A(eqn, opts, eqn.type, z, 'N');
y2 = oper.mul_A(eqn, opts, adjoint, x, 'N');

%% Bilinear Term
numberOf_N_matrices = length(eqn.N_);
rowN = size(eqn.N_{1}, 1);
colZ = size(Z, 2);
bilinearSUM = zeros(rowN, numberOf_N_matrices * colZ);

colCompress_start = 1;
colCompress_end = colZ;

% Building [N_1*Z, ..., N_k*Z] with N_k' for type 'T'
for currentN_k = 1:numberOf_N_matrices

    bilinearSUM(:, colCompress_start:colCompress_end) = ...
        oper.mul_N(eqn, opts, eqn.type, Z, 'N', currentN_k);

    colCompress_start = colCompress_end + 1;
    colCompress_end = colCompress_end + colZ;
end

switch eqn.type
    case 'N'
        y3 = bilinearSUM * (bilinearSUM' * x);
    case 'T'
        y3 = bilinearSUM * (bilinearSUM' * x);
    otherwise
        mess_err(opts, 'control_data', 'eqn.type has to be ''N'' or ''T''');
end

%% Terms with A
y2 = Z * (Z' * y2);

if eqn.haveE
    y2 = oper.mul_E(eqn, opts, eqn.type, y2, 'N');
end

%% Complete SUM
y = y1 + y2 + y3;

switch eqn.type
    case 'N'
        y = y + eqn.pB * (eqn.pB' * x);
    case 'T'
        y = y + eqn.pC' * (eqn.pC * x);
    otherwise
        mess_err(opts, 'control_data', 'eqn.type has to be ''N'' or ''T''');
end

end
