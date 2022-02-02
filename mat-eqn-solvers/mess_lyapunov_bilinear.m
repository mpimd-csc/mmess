function [out, eqn, opts, oper] = mess_lyapunov_bilinear(eqn, opts, oper)
% Solve continuous-time Lyapunov plus positive equations with sparse
% coefficients
%
%         A*P*E' + E*P*A' + Sum_N_k*P*N_k' + B*B' = 0       (N)
%         A'*Q*E + E'*Q*A + Sum_N_k'*Q*N_k + C'*C = 0       (T)
%
% Algorithm gives you Z such that Z*Z' approximates P, or Q respectively.
%
% Currently only 'default' usfs are supported.
%
% Input
%   eqn                 struct containing data for equations
%
%   opts                struct containing parameters for the algorithm and
%                       subalgorithms (see below)
%
%   oper                struct contains function handles for operation
%                       with A, E and N
%
% Output
%   out                 struct containing output information, especially
%                       out.Z is the solution factor
%
%   eqn, opts, oper     same as inputs
%
% Input fields in struct eqn:
%   eqn.B       dense (n x m1) matrix B
%
%   eqn.C       dense (m2 x n) matrix C
%
%   eqn.N_      for the 'default' usfs: Cell array with N_k = N{k} for
%               k = 1,2, ...
%               (if all N_k are given in one large matrix  of size
%                n x (n * number_of_matrices)
%               it will be transformed to a cell and restored on the
%               output)
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) above is solved
%               (optional, default fallback: 'N')
%
%   eqn.haveE   possible  values: 0, 1, false, true
%               if haveE = 0: matrix E is assumed to be the identity
%               (optional, default: 0)
%
%   Depending on the operator chosen by the operatormanager, additional
%   fields may be needed. For the "default", e.g., eqn.A_ and eqn.E_ hold
%   the A and E matrices.
%
%
% Input fields in struct opts:
%   opts.blyap.res_tol          possible  values: real > 0
%                               residual norm tolerance
%                               (optional, default: 1e-5)
%
%   opts.blyap.rel_diff_tol     possible  values: real > 0
%                               stagnation tolerance for the relative
%                               difference of subsequent iterates
%                               (optional, default: 1e-5)
%
%   opts.blyap.maxiter          possible  values: integer > 0
%                               maximum iteration number
%                               (optional, default: 10)
%
%   opts.blyap.info             possible values 0, 1
%                               dis-/enables additional output information in
%                               each iteration (current step, residual, rank,
%                               elapsed time)
%                               (optional, default: 0)
%
%   opts.adi                    options structure can be used to pass setting
%                               to the LRADI or ADI shift computation(optional)
%                               (see corresponding routines for additional
%                                information)
%
%   opts.resopts                options structure can be used to pass setting
%                               to the mess_res2_norm computation (optional)
%                               (see corresponding routines for additional
%                                information)
%
%
% Input fields in struct opts:
%   oper                        = operatormanager('default')
%                               (other cases not implemented yet)
%
%
% Output fields in struct out:
%   out.Z                       with out.Z*outZ' = P - Error
%
%
%   out.niter_bilinear          number of mess_lyapunov_bilinear iterations
%
%
%   out.resNorm_bilinear        Residual Value as in
%                               ||A*P'*E' + E*P*A' + Sum_N_k*P*N_k' + B*B||
%                                with P=Z*Z'
%                               (NOT calculated in every iteration, but
%                               only when relative change is small enough,
%                               or rank of Z stagnates; 0 when not computed)
%
%   out.rankZ_bilinear          array containing rank(Z) of every Iteration
%
% References
%
% [1] P. K. Goyal, System-theoretic model order reduction for bilinear and
%     quadratic-bilinear systems, Dissertation, Department of Mathematics,
%     Otto von Guericke University, Magdeburg, Germany (2018).
%     https://doi.org/10.25673/5319.
%
% [2] P. Benner, P. Goyal, Balanced truncation model order reduction for
%     quadratic-bilinear systems, e-prints 1705.00160, arXiv, math.OC (2017).
%     URL https://arxiv.org/abs/1705.00160
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


%% Check oper
if not(isequal(oper.name, 'default'))
    error('MESS:notimplemented', ...
          ['Feature not yet implemented!', ...
           'Only accepts operatormanager(''default'')']);
end

%% Check for eqn
if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
if not(isfield(eqn, 'type')), eqn.type = 'N'; end

% Check for A_ and E_
[~, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');

% Check for B and C
if not(isfield(eqn, 'B')) || isempty(eqn.B)
    warning('MESS:control_data','eqn.B is missing');
end

if not(isfield(eqn, 'C')) || isempty(eqn.C)
    warning('MESS:control_data','eqn.C is missing');
end

% Check for eqn.type and declare B or C as full
switch eqn.type
    case 'N'
        eqn.B = full(eqn.B);
    case 'T'
        eqn.C = full(eqn.C);
    otherwise
        error('MESS:control_data','eqn.type has to be ''N'' or ''T''');
end

% Transform N to Cell if given as matrix and check N
[eqn, opts, oper] = oper.mul_N_pre(eqn, opts, oper);

%% check opts for substructures
if not(isfield(opts, 'norm')), opts.norm ='fro'; end

expected_subStructures = {'blyap', 'shifts','adi', 'fopts'};

expected_sub_subStructures{1} = {'maxiter', 'rel_diff_tol', 'res_tol'};
default_values{1} = {10, 1e-5, 1e-5};

expected_sub_subStructures{2} = {'num_Ritz', 'num_hRitz', ...
                                 'method', 'num_desired'};
default_values{2} = {50, 25, 'projection', 6};

expected_sub_subStructures{3} = {'maxiter', 'res_tol', ...
    'rel_diff_tol', 'norm'};
default_values{3} = {100, 1e-12, 0, 'fro'};

expected_sub_subStructures{4} = {'LDL_T', 'norm'};
default_values{4} = {0, 'fro'};

opts = check_opts(opts, expected_subStructures, expected_sub_subStructures, ...
                  default_values);

% check res2_norm opts
if not(isfield(opts, 'resopts')) || isempty(opts.resopts)
    warning('MESS:control_data',...
            'opts.resopts is missing. parameters set to default values');
    fields = {'res'};
    c = cell(length(fields), 1);
    opts.resopts = cell2struct(c, fields);
end

expected_subStructures_resopts = {'res'};
expected_sub_subStructures_resopts{1} = {'maxiter', 'tol', 'orth'};
default_values_resopts{1} = {10, 1e-6, 0};

opts.resopts = check_opts(opts.resopts,...
    expected_subStructures_resopts,expected_sub_subStructures_resopts,...
    default_values_resopts);

%% Initialize Data
eqn_iter = eqn; % Equationset per Iteration
B_iter = eqn.B; % changes per Iteration
C_iter = eqn.C; % changes per Iteration

Z = [];
rz_new = 0;

ranksZ = zeros(opts.blyap.maxiter, 1); % to save output
lyapNorm = zeros(opts.blyap.maxiter, 1); % to save output

rowA = size(eqn.A_, 1);
numberOf_N_matrices = length(eqn.N_);

switch eqn.type
    case 'N'
        norm_B_or_C = norm((eqn.B' * eqn.B));
    case 'T'
        norm_B_or_C = norm(eqn.C * eqn.C');
    otherwise
        error('MESS:control_data','eqn.type has to be ''N'' or ''T''');
end

% info opts
if not(isfield(opts.blyap, 'info'))
    opts.blyap.info = 0;
else
    if not(isnumeric(opts.blyap.info)) && not(islogical(opts.blyap.info))
        error('MESS:info',...
              'opts.shifts.info parameter must be logical or numeric.');
    end
    tic
end

%% Start iteration
for k = 1 : (opts.blyap.maxiter)
    % Solve A*VV'E' + EVV'A'  + BB' = 0 or
    % Solve A'VV'E  + E'VV'A  + C'C = 0
    eqn_iter.B = B_iter;
    eqn_iter.C = C_iter;

    out_lradi = mess_lradi(eqn_iter, opts, oper);
    V =  mess_column_compression( out_lradi.Z, 'N', [], eps);

    % column_compression
    colV = size(V, 2);
    compress = zeros(rowA, numberOf_N_matrices * colV);

    col_start = 1;
    % building [N_1*V, ..., N_k*V], or the same with N_k' for type 'T'
    for currentN_k = 1:numberOf_N_matrices
        NV = mess_column_compression(oper.mul_N(eqn, opts, eqn.type, V, 'N', ...
                                     currentN_k), 'N', [], eps);
        col_NV = size(NV, 2);
        columns = col_start:(col_start+col_NV-1);
        compress(:, columns) = NV;
        col_start = col_start + col_NV;
    end

    switch eqn.type
        case 'N'
            B_iter = mess_column_compression(compress(:, 1:col_start-1), ...
                                             'N', [], eps);
        case 'T'
            C_iter = mess_column_compression(compress(:, 1:col_start-1), ...
                                             'N', [], eps);
            C_iter = C_iter';
        otherwise
            error('MESS:control_data','eqn.type has to be ''N'' or ''T''');
    end

    Z = mess_column_compression ([Z,V], 'N', [], eps);

    % calculate exit conditions
    rZ_old = rz_new;
    rz_new = rank(Z);

    ranksZ(k) = rank(Z);


    normZ = norm(Z' * Z);
    normV = norm(V)^2;

    % test the conditions
    if ((normV/normZ) < opts.blyap.rel_diff_tol) || (rZ_old == rz_new)
        % calculate ||A*Z*Z'*E' + E*Z*Z'*A' + Sum_N_k*P*N_k' + B*B|| for
        % current Z
        [lyapunov_Norm, ~, ~, ~, eqn, opts.fopts, oper] = mess_res2_norms(Z,...
            'lyapunov_QB', eqn, opts.fopts, oper, opts.resopts);

        lyapNorm(k) = lyapunov_Norm;

        if lyapunov_Norm < (opts.blyap.res_tol * norm_B_or_C)
            niter_bilinear = k;

            if opts.blyap.info
                fprintf('step: %4d residual: %e rank(Z): %d \n', ...
                        k, lyapNorm(k), rz_new);
                fprintf(['Converged at step: %4d with residual: %e rank(Z):', ...
                        '%d \n'], k, lyapNorm(k), rz_new);
            end
            break
        end
    end

    if k == opts.blyap.maxiter
        [lyapunov_Norm, ~, ~, ~, eqn, opts.fopts, oper] = ...
            mess_res2_norms(Z,'lyapunov_QB', eqn, opts.fopts, ...
            oper, opts.resopts);

        lyapNorm(k) = lyapunov_Norm;
        niter_bilinear = k;

        warning('MESS:maxiter', ...
                'REACHED MAXIMUM ITERATION NUMBER, FINAL RESIDUAL: %d', ...
                lyapunov_Norm);

        if opts.blyap.info
            fprintf('step: %4d residual: %e rank(Z): %d \n', ...
                    k, rz_new, lyapNorm(k));
            fprintf(['Did not converge after %4d steps last residual:', ...
                     ' %e rank(Z): %d \n'], k, lyapNorm(k), rz_new);
        end
    end

    if opts.blyap.info
        fprintf('step: %4d residual: %e rank(Z): %d \n', ...
                k, lyapNorm(k), rz_new);
    end
end

% Elapsed time
if opts.blyap.info
    toc
end

% Saving output
out.Z = Z;

out.niter_bilinear = niter_bilinear;

% format rank output
out.rankZ_bilinear = zeros(out.niter_bilinear, 1);
for k = 1 : out.niter_bilinear
    out.rankZ_bilinear(k) = ranksZ(k);
end

% format residual output
out.resNorm_bilinear = zeros(out.niter_bilinear, 1);
for k = 1 : out.niter_bilinear
    out.resNorm_bilinear(k) = lyapNorm(k);
end

% Transform N back into Matrix (only if it wasn't given as a cell)
[eqn, opts, oper] = oper.mul_N_post(eqn, opts, oper);

end

%% Local function: check_opts
function opts = check_opts(opts, expected_subStructures, ...
                           expected_sub_subStructures, default_values)
% helper function

% Check for opts
% Input
%   opts                             struct to be checked for options
%
%   expected_subStructures           input as cell with the names of the
%                                    necessary substructures
%
%   expected_sub_subStructures       input as cell with the names of the
%                                    necessary sub_substructures
%                                    {n} names of the nth substructure
%
%   default_values                   input as cell with the default values for
%                                    the necessary sub_substructures
%                                    {n} values of the nth substructure
%
% Output
%   opts                             opts structure with the new substructures
%                                    using default values if not assigned
%
%

for k = 1 : length(expected_subStructures)
    % checks if all necessary data is given
    if not(length(expected_subStructures) == length(expected_sub_subStructures))
        error('MESS:control_data', ...
              'Sub_SubStructure for some SubStructure is missing!')
    end

    if not(length(expected_sub_subStructures{k}) == length(default_values{k}))
        error('MESS:control_data', ...
              'Not every sub_subStructure has an assigned default value!');
    end

    % creates structures if not yet existing
    if not(isfield(opts, (expected_subStructures{k}))) || ...
       isempty(opts.(expected_subStructures{k}))
        warning('MESS:control_data', ...
                'opts.%s is missing. parameters set to default values', ...
                expected_subStructures{1});

        create_empty_opt = cell(length(expected_sub_subStructures{k}), 1);
        opts.(expected_subStructures{k}) = cell2struct(create_empty_opt, ...
                                               expected_sub_subStructures{k});
    end

    for l = 1 : length(expected_sub_subStructures{k})

        % assigins default values if there are no given values
        if not(isfield(opts.(expected_subStructures{k}), ...
               expected_sub_subStructures{k}{l})) || ...
           isempty(opts.(expected_subStructures{k}).(expected_sub_subStructures{k}{l}))

            default_message = default_values{k}{l};

            if not(ischar(default_message))
                default_message = num2str(default_values{k}{l});
            end

            warning('MESS:control_data','opts.%s.%s is set to %s (default)', ...
                    expected_subStructures{k}, expected_sub_subStructures{k}{l}, default_message);

            opts.(expected_subStructures{k}).(expected_sub_subStructures{k}{l}) = default_values{k}{l};
        end
    end
end

end
