function [IQLas, IQDas, eqn, opts, oper] ...
    = IQ(eqn, opts, oper, t0s, h, as)
% Compute the low-rank factors LDL^T of the integrals
%
% \int_{t0s(i)}^{t0s(i) + as(k)*h}{
%     exp(s(A E^{-1})^T)E^{-T} C^T C E^{-1} exp(s A E^{-1}) ds},
%
% for all the initial times t0(i) in the vector t0s and all scaling factors
% as(k) in the vector as. This is for eqn.type == 'T'. For eqn.type == 'N',
% we factorize instead
%
% \int_{t0s(i)}^{t0s(i) + h}{
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
%   t0s     starting times, array
% 
%   h       main time step size, scalar > 0
%
%   as      list of a factors, array
%
%
% Output:
%   IQLas   cell array such that IQLAs{i}{k} is the L factor corresponding
%           to t0s(i) and as(k)
%
%   IQDas   cell array such that IQLAs{i}{k} is the D factor corresponding
%           to t0s(i) and as(k)
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009-2019
%

type = opts.splitting.quadrature.type;
order = opts.splitting.quadrature.order;

% Decide which type of quadrature to use and compute nodes and weights for
% the interval [0, h]. These are then shifted appropriately later.
adaptive = false;
if strcmp(type, 'gauss')
    [xj, wj] = gauss_quadrature_parameters(h, order);
elseif strcmp(type, 'clenshawcurtis')
    if rem(order,2) == 1
        error('MESS:control_data', ['Quadrature order must be even ' ...
                            'Clenshaw-Curtis.']);
    end
    [xj, wj] = clenshawcurtis_parameters(0, h, order);
elseif strcmp(type, 'adaptive')
    % nodes/weights computed in loop further down
    adaptive = true;
    order = 4;
elseif strcmp(type, 'equidistant')
    % Equidistant nodes - used for the adaptive methods
    xj = linspace(0, h, order+1)';
    wj = compute_quadrature_weights(h, xj);
end

TOL = 1e-4;
if isfield(opts.splitting.quadrature, 'tol')
    TOL = opts.splitting.quadrature.tol;
end

nts = length(t0s);
nas = length(as);
IQLas = cell(nts, nas);
IQDas = cell(nts, nas);
if adaptive
    IQLas2 = cell(nts, nas);
    IQDas2 = cell(nts, nas);
end
Las = cell(nts,nas); % The Ls used to form each IQL

% In the autonomous case, we have a fixed block to apply the matrix
% exponential to:
if not(eqn.LTV)
    if eqn.type == 'T'
        LQ = eqn.C';
    elseif eqn.type == 'N'
        LQ = eqn.B;
    end
    % We assume that DQ = I, always.
    DQ = eye(size(LQ, 2));
end


for i = 1:length(t0s)
    t0 = t0s(i);
    errest = Inf;
    % This loop breaks after one iteration unless the adaptive strategy is
    % used:
    while errest > TOL 
        if adaptive
            [xj, wj] = clenshawcurtis_parameters(0, h, order);
            [~, wj2] = clenshawcurtis_parameters(0, h, order/2);
            % Sort the nodes from 0 to h instead, for iterative computation
            % (the weights are symmetric)
            xj = flipud(xj);
        end
        
        % For each scaling factor as(k) compute the values Ls(xj)
        for k = 1:length(as)
            if not(eqn.LTV) % Autonomous case
                Ls = cell(1, length(xj));
                RHS = oper.sol_E(eqn, opts, eqn.type, LQ, 'N');
                % For legacy reasons, the sol_E_DAE2 function does not 
                % return the rows corresponding to only the states, but the 
                % full solution. We therefore cut it down to size here.
                % This does nothing in the other cases when RHS already has
                % the correct size.
                RHS = RHS(1:size(LQ,1), :); % 
                [out, ~, opts, ~] ...
                    = mess_exp_action(eqn, opts, oper, as(k)*xj(1), RHS);
                Ls{1} = out.Z;
                for j = 2:length(xj)
                    dx = as(k)*(xj(j) - xj(j-1));
                    Ltemp = oper.mul_E(eqn, opts, eqn.type, Ls{j-1}, 'N');
                    RHS = oper.sol_E(eqn, opts, eqn.type, Ltemp, 'N');
                    % See above comment about RHS:
                    RHS = RHS(1:size(LQ,1), :); 
                    [out, ~, opts, ~] ...
                        = mess_exp_action(eqn, opts, oper, dx, RHS);
                    Ls{j} = out.Z;
                end
            else % Timevarying case
                % This can be done somewhat faster by using BLAS-3 
                % blocking, but it requires more memory. See the DREsplit 
                % package for this option, omitted in the M.E.S.S. version.
                for j = 1:length(xj)
                    t = t0 + as(k)*xj(j);
                    % RHS is E(t)'\LQ(t), so update the matrices
                    [eqn, opts, oper] = ...
                        opts.splitting.eval_matrix_functions(eqn, opts, ...
                                                             oper, t);
                    if eqn.type == 'T'
                        LQ = eqn.C';
                    elseif eqn.type == 'N'
                        LQ = eqn.B;
                    end
                    ELQ = oper.sol_E(eqn, opts, eqn.type, LQ, 'N');
                    % See above comment about RHS:
                    ELQ = ELQ(1:size(LQ,1), :); 
                    
                    
                    [out, ~, opts, ~] ...
                        = mess_exp_action(...
                            eqn, opts, oper, ...
                            as(k)*(h-xj(j)), ...
                            full(ELQ), ...
                            t);
                    Ls{j} = out.Z;
                end
            end
            Las{i}{k} = Ls;
            IQLas{i}{k} = cell2mat(Ls);
            % Build block diagonal matrix with weight(j)*as(k)*DQ as blocks
            if eqn.LTV % We don't know the size of DQ a priori.
                DQ = eye(size(LQ, 2));
            end
            IQDas{i}{k} = kron(diag(wj * as(k)), DQ); 
        end
        % Column compress approximations
        for k = 1:length(as)
            [IQLas{i}{k}, IQDas{i}{k}] = ...
                mess_column_compression(IQLas{i}{k}, 'N', IQDas{i}{k}, ...
                opts.splitting.trunc_tol, opts.splitting.trunc_info);
        end
        
        if adaptive
            % Column compress coarser approximation
            for k = 1:length(as)
                [IQLas2{i}{k}, IQDas2{i}{k}] = ...
                    mess_column_compression(...
                        cell2mat(Las{i}{k}(1:2:end)), 'N', ...
                        kron(diag(wj2 * as(k)), DQ), ...
                        opts.splitting.trunc_tol, ...
                        opts.splitting.trunc_info);
            end
            errest = 0;
            for k = 1:length(as)
                errest = ...
                max(errest, ...
                    outerfrobnormdiff_LDLT(IQLas{i}{k}, IQDas{i}{k}, ...
                                           IQLas2{i}{k}, IQDas2{i}{k}));
            end
        else
            break
        end
        order = order * 2;
    end
end

end