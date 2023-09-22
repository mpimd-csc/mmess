function [IQLas, IQDas, eqn, opts, oper] ...
    = IQ(eqn, opts, oper, t_zeros, h, as)
% Compute the low-rank factors LDL^T of the integrals
%
% \int_{t_zeros(i)}^{t_zeros(i) + as(k)*h}{
%     exp(s(A E^{-1})^T)E^{-T} C^T C E^{-1} exp(s A E^{-1}) ds},
%
% for all the initial times t0(i) in the vector t_zeros and all scaling factors
% as(k) in the vector as. This is for eqn.type == 'T'. For eqn.type == 'N',
% we factorize instead
%
% \int_{t_zeros(i)}^{t_zeros(i) + h}{
%     exp(s A E^{-1})E^{-1} B B^T E^{-T} exp(s (A E^{-1})^T) ds}.
%
% If eqn.LTV is true, the integrand is replaced by the relevant integral.
%
%
% Input & Output
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   oper            struct contains function handles for operations
%                   with A and E
%
% Input:
%   t_zeros starting times, array
%
%   h       main time step size, scalar > 0
%
%   as      list of a factors, array
%
%
% Output:
%   IQLas   cell array such that IQLAs{i}{k} is the L factor corresponding
%           to t_zeros(i) and as(k)
%
%   IQDas   cell array such that IQLAs{i}{k} is the D factor corresponding
%           to t_zeros(i) and as(k)
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

qtype = opts.splitting.quadrature.type;
order = opts.splitting.quadrature.order;
rel = opts.splitting.quadrature.rel_err;

% Decide which type of quadrature to use and compute nodes and weights for
% the interval [0, h]. These are then shifted appropriately later.
adaptive = opts.splitting.quadrature.adaptive;
embedded = adaptive && opts.splitting.quadrature.embedded;
intervals = opts.splitting.quadrature.intervals;

if (strcmp(qtype, 'clenshawcurtis') || embedded) && rem(order, 2) == 1
    mess_err(opts, 'control_data', ['Quadrature order must be even for' ...
                                    'Clenshaw-Curtis or embedded formulas.']);
end
switch qtype
    case 'gauss'
        [xj_base, wj_base] = gauss_quadrature_parameters(1, order);
        embedded = false;
    case 'clenshawcurtis'
        [xj_base, wj_base] = clenshawcurtis_parameters(0, 1, order);
        xj_base = flipud(xj_base);
        if embedded
            [~, wj2_base] = clenshawcurtis_parameters(0, 1, order / 2);
        end
    case 'equidistant'
        xj_base = linspace(0, 1, order + 1)';
        wj_base = compute_quadrature_weights(1, xj_base);
        if embedded
            wj2_base = compute_quadrature_weights(1, xj_base(1:2:end));
        end
    otherwise
        mess_err(opts, 'control_data', ...
                 'The specified quadrature type is not supported.');
end

%% set data
TOL = 1e-4;

if isfield(opts.splitting.quadrature, 'tol')
    TOL = opts.splitting.quadrature.tol;
end

nts = length(t_zeros);
nas = length(as);
IQLas = cell(nts, nas);
IQDas = cell(nts, nas);

% In the autonomous case, we have a fixed block to apply the matrix
% exponential to:
if not(eqn.LTV)

    if eqn.type == 'T'
        LQ = eqn.C';
    elseif eqn.type == 'N'
        LQ = eqn.B;
    end
else
    LQ = [];
end
%% iterate t_zeros
for i = 1:nts
    t0 = t_zeros(i);
    %% For each scaling factor as(k) approximate the integral over
    % [t0, t0 + as(k)*h]
    % Start at the largest a, then the final number of intervals should
    % be enough also for the smaller ones.
    for k = 1:nas
        errest = Inf;
        IQLf = [];
        %% adaptive loop
        % This loop breaks after one iteration unless the adaptive strategy
        % is used:
        while errest > TOL
            % Coarse intervals
            % Scale from [0,1] to [0, as(k)*h / intervals]
            interval_length = as(k) * h / intervals;
            xjc = xj_base * interval_length;
            wjc = wj_base * interval_length;

            % The intervals are t0 + xjc, t0 + interval_length + xjc, ...
            % Starting points of these intervals:
            t0c = t0 + interval_length * (0:intervals); % one extra

            % Fine intervals, twice as many
            interval_length_f = as(k) * h / intervals / 2;
            xjf = xj_base * interval_length_f;
            t0f = t0 + interval_length_f * (0:intervals * 2); % one extra
            wjf = wj_base * interval_length_f;
            if embedded
                wj2f = wj2_base * interval_length_f;
            else
                wj2f = [];
            end
            %% build block matrices L and D
            if isempty(IQLf)
                [eqn, opts, oper, IQLc, IQDc] = ...
                    approximate_subintegrals(eqn, opts, oper, t0c, xjc, ...
                                             wjc, [], as(k), LQ, false, rel);
            else
                IQLc = IQLf;
                IQDc = IQDf;
            end
            [eqn, opts, oper, IQLf, IQDf, errest_embedded] = ...
                approximate_subintegrals(eqn, opts, oper, t0f, xjf, wjf, ...
                                         wj2f, as(k), LQ, embedded, rel);

            % Final result (temporary)
            IQLas{i}{k} = IQLf;
            IQDas{i}{k} = IQDf;

            %% error estimation
            errest_composite = outerfrobnormdiff_LDLT( ...
                                                      IQLc, IQDc, IQLf, IQDf, rel);
            errest = max(errest_embedded, errest_composite);
            if adaptive
                if errest > TOL
                    intervals = intervals * 2;
                end
            else
                break
            end
        end
    end
end
end

function [eqn, opts, oper, IQL, IQD, errest] = ...
    approximate_subintegrals(eqn, opts, oper, t_zeros_subintegrals, ...
                             xj, wj, wj2, a, LQ, embedded, rel)
IQL = cell(1, length(t_zeros_subintegrals) - 1);
IQD = cell(1, length(t_zeros_subintegrals) - 1);
IQL2 = cell(1, length(t_zeros_subintegrals) - 1);
IQD2 = cell(1, length(t_zeros_subintegrals) - 1);
Ls = cell(1, length(xj));
if not(eqn.LTV)
    if eqn.type == 'T'
        LQ = eqn.C';
    elseif eqn.type == 'N'
        LQ = eqn.B;
    end
end
for l = 1:length(t_zeros_subintegrals) - 1
    for j = 1:length(xj)
        r = t_zeros_subintegrals(l) + xj(j);
        t = t_zeros_subintegrals(end);
        if eqn.LTV
            % RHS is E(r)'\LQ(r), so update the matrices
            [eqn, opts, oper] = ...
                opts.splitting.eval_matrix_functions(eqn, opts, ...
                                                     oper, r);
            if eqn.type == 'T'
                LQ = eqn.C';
            elseif eqn.type == 'N'
                LQ = eqn.B;
            end
        end
        ELQ = oper.sol_E(eqn, opts, eqn.type, LQ, 'N');
        % See above comment about RHS:
        ELQ = ELQ(1:size(LQ, 1), :);

        [out, ~, opts] = mess_exp_action(eqn, opts, oper, ...
                                         t - r, full(ELQ), r);
        Ls{j} = out.Z;
    end

    %% build block matrices L and D over subinterval
    IQL{l} = cell2mat(Ls);
    % Build block diagonal matrix with wj(j)*DQ as blocks
    % LTV: We don't know the size of DQ a priori.
    DQ = eye(size(LQ, 2));
    IQD{l} = kron(diag(wj), DQ);

    % Column compress approximation
    [IQL{l}, IQD{l}] = ...
        mess_column_compression(IQL{l}, 'N', IQD{l}, ...
                                opts.splitting.trunc_tol, opts.splitting.trunc_info);
    % Column compress coarser approximation
    if embedded
        IQL2{l} = cell2mat(Ls(1:2:end));
        IQD2{l} = kron(diag(wj2 * a), DQ);
        [IQL2{l}, IQD2{l}] = ...
            mess_column_compression(IQL2{l}, 'N', IQD2{l}, ...
                                    opts.splitting.trunc_tol, opts.splitting.trunc_info);
    end
end

%% build block matrices L and D over full interval
% Column compress over t0 to get approximation over full interval
[IQL, IQD] = ...
    mess_column_compression(cell2mat(IQL), 'N', blkdiag(IQD{:}), ...
                            opts.splitting.trunc_tol, opts.splitting.trunc_info);
if embedded
    % Column compress coarser approximation
    [IQL2, IQD2] = ...
        mess_column_compression(cell2mat(IQL2), 'N', blkdiag(IQD2{:}), ...
                                opts.splitting.trunc_tol, opts.splitting.trunc_info);
    %% error estimation
    errest = outerfrobnormdiff_LDLT(IQL, IQD, IQL2, IQD2, rel);
else
    errest = 0;
end
end
