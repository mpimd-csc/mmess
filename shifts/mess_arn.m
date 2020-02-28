function [H,V] = mess_arn(eqn, opts, oper, opA)
%
%  Arnoldi method w.r.t. opA(A)
%
%  Calling sequence:
%
%    [H,V] = mess_arn(eqn, opts, oper, opA)
%
%  Input:
%
%    eqn       eqn contains data for A (Input matrix) for which
%              the Arnoldi Algorithm is supposed to run.
%    opts      struct contains parameters for the algorithm
%    oper      contains function handles with operations for A and E
%    opA       character specifies the form of opA(A)
%              opA = 'N' performs Arnoldi method
%              opA = 'I' performs inverse Arnoldi method
%
%  Output:
%
%    H         matrix H ((k+1)-x-k matrix, upper Hessenberg);
%    V         matrix V (n-x-(k+1) matrix, orthogonal columns).
%
%  Method:
%
%    The Arnoldi method produces matrices V and H such that
%
%      V(:,1) in span{b0},
%      V'*V = eye(k+1),
%      F*V(:,1:k) = V*H.
%
%    b0 = opts.shifts.b0
%
%  Remark:
%
%    This implementation does not check for (near-)breakdown!
%
%   uses operatorfunctions size, sol_A, mul_A, sol_E, mul_E

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
%               2009-2020
%


% Input data not completely checked!e
%% check input data
if not(isfield(opts,'shifts')) || not(isstruct(opts.shifts))
    warning('MESS:control_data',['shift parameter control structure missing.', ...
        'Switching to default num_desired = 25, num_Ritz = 50, num_hRitz = 25.']);
    opts.shifts.num_desired = 25;
    opts.shifts.num_Ritz = 50;
    opts.shifts.num_hRitz = 25;
else
    if not(isfield(opts.shifts,'num_desired'))||not(isnumeric(opts.shifts.num_desired))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_desired field.', ...
            'Switching to default: 25']);
        opts.shifts.num_desired = 25;
    end
    if strcmp(opts.shifts.method,'heur')&&...
       (not(isfield(opts.shifts,'num_Ritz'))||not(isnumeric(opts.shifts.num_Ritz)))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_Ritz field.', ...
            'Switching to default: 50']);
        opts.shifts.num_Ritz = 50;
    end
    if strcmp(opts.shifts.method,'heur')&&...
       (not(isfield(opts.shifts,'num_hRitz'))||not(isnumeric(opts.shifts.num_hRitz)))
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.shifts.num_hRitz field.', ...
            'Switching to default: 25']);
        opts.shifts.num_hRitz = 25;
    end
end

if not(isfield(eqn, 'haveE')), eqn.haveE = 0; end
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');
if not(result)
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end
if isfield(eqn, 'haveUV') && eqn.haveUV
    if not(isfield(eqn,'U')) || isempty(eqn.U) || not(isfield(eqn,'V')) || isempty(eqn.V)...
            || not((size(eqn.U,1))==size(eqn.V,1) && size(eqn.U,2)==size(eqn.V,2))
        error('MESS:SMW','Inappropriate data of low rank updated opertor (eqn.U and eqn.V)');
    end
end
if not(isfield(opts,'rosenbrock')), opts.rosenbrock=[]; end
if isstruct(opts.rosenbrock)&&isfield(opts.rosenbrock,'tau')
    rosenbrock = 1;
    pc = -1 / (2 * opts.rosenbrock.tau);
else
    rosenbrock = 0;
end
if not(isfield(opts,'bdf')), opts.bdf=[]; end
if isstruct(opts.bdf) && isfield(opts.bdf, 'tau') && isfield(opts.bdf, 'beta')
    bdf = 1;
    pc = -1 / (2 * opts.bdf.tau * opts.bdf.beta);
else
    bdf = 0;
end

%% check input Paramters
if not(ischar(opA))
    error('MESS:error_arguments', 'opA is not a char');
end

opA = upper(opA);
if(not((opA == 'N' || opA == 'I')))
    error('MESS:error_arguments','opA is not ''N'' or ''I''');
end

% returns order of A or states of A, A is supposed to be square
n = oper.size(eqn, opts);

if opA == 'I'
    k = opts.shifts.num_hRitz;
else
    k = opts.shifts.num_Ritz;
end

if k >= n - 1, error('k must be smaller than the order of A!'); end
%% initialize data
if (not(isfield(opts.shifts, 'b0')) || isempty(opts.shifts.b0))
    b0 = ones(n,1);
else
    b0 = opts.shifts.b0;
end

V = zeros(length(b0), k + 1);
H = zeros(k + 1, k);

V(:, 1) = (1.0 / norm(b0)) * b0;

beta = 0;
%% perform Arnoldi method
for j = 1 : k
    
    if j > 1
        V(:, j) = (1.0 / beta) * w;
    end
       
    % no eqn.type cases needed, eigenvalues are the same for transposed
    % operator
    if opA == 'I' %Perform inverse Arnodi
        if isfield(eqn,'haveE') && eqn.haveE
            if isfield(eqn,'haveUV') && eqn.haveUV
                if bdf
                    AB = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                        1 / (opts.bdf.tau * opts.bdf.beta) * ...
                        [oper.mul_E(eqn, opts, 'N', V(:, j), 'N') eqn.U], 'N');
                elseif rosenbrock
                    AB = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                        [oper.mul_E(eqn, opts, 'N', V(:, j), 'N') eqn.U], 'N');
                else
                    AB = oper.sol_A(eqn, opts, 'N', [oper.mul_E(eqn, opts, 'N', V(:, j), 'N') eqn.U], 'N');
                end
                AV = AB(:,1);
                AU = AB(:, 2 : end);
                w = AV - AU * ((speye(size(eqn.U, 2)) + eqn.V' * AU) \ (eqn.V' * AV));
            else
                if bdf
                    w = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                        1 / (opts.bdf.tau * opts.bdf.beta) * ...
                        oper.mul_E(eqn, opts, 'N', V(:, j), 'N'), 'N');
                else
                    w = oper.sol_A(eqn, opts, 'N', oper.mul_E(eqn, opts, 'N', V(:, j), 'N'), 'N');
                end
            end
        else
            if isfield(eqn,'haveUV') && eqn.haveUV
                if bdf
                    AB = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                        1 / (opts.bdf.tau * opts.bdf.beta) * ...
                        [V(:, j) eqn.U], 'N');
                elseif rosenbrock
                    AB = oper.sol_ApE(eqn, opts, 'N', pc, 'N', [V(:, j) eqn.U], 'N');
                else
                    AB = oper.sol_A(eqn, opts, 'N', [V(:, j) eqn.U], 'N');
                end
                AV = AB(:,1);
                AU = AB(:, 2 : end);
                w = AV - AU * ((speye(size(eqn.U, 2)) + eqn.V' * AU) \ (eqn.V' * AV));
            else
                if bdf
                    w = oper.sol_ApE(eqn, opts, 'N', pc, 'N', ...
                        1 / (opts.bdf.tau * opts.bdf.beta) * ...
                        V(:, j), 'N');
                else
                    w = oper.sol_A(eqn, opts, 'N', V(:, j), 'N');
                end
            end
        end
    else % opA = 'N' Perform standard Arnoldi
        if isfield(eqn,'haveE') && eqn.haveE
            if isfield(eqn,'haveUV') && eqn.haveUV
                if bdf
                    w = oper.sol_E(eqn, opts, 'N',...
                        oper.mul_ApE(eqn, opts, 'N', pc, 'N', ...
                        (opts.bdf.tau * opts.bdf.beta) * ...
                        V(:, j), 'N') + eqn.U * (eqn.V' * V(:, j)), 'N');
                elseif rosenbrock
                    w = oper.sol_E(eqn, opts, 'N',...
                        oper.mul_ApE(eqn, opts, 'N', pc, 'N', V(:, j), 'N')...
                        + eqn.U * (eqn.V' * V(:, j)), 'N');
                else
                    w = oper.sol_E(eqn, opts, 'N',...
                        oper.mul_A(eqn, opts, 'N', V(:, j), 'N')...
                        + eqn.U * (eqn.V' * V(:, j)), 'N');
                end
            else
                if bdf
                    w = oper.sol_E(eqn, opts, 'N', (opts.bdf.tau * opts.bdf.beta) * ...
                        oper.mul_ApE(eqn, opts, ...
                        'N', pc, 'N', V(:, j), 'N'), 'N');
                else
                    w = oper.sol_E(eqn, opts, 'N', oper.mul_A(eqn, opts, 'N', V(:, j), 'N'), 'N');
                end
            end
        else
            if isfield(eqn,'haveUV') && eqn.haveUV
                if bdf
                    w = (opts.bdf.tau * opts.bdf.beta) * ...
                        oper.mul_ApE(eqn, opts, 'N', pc, 'N', V(:, j), 'N') ...
                        + eqn.U * (eqn.V' * V(:, j));
                elseif rosenbrock
                    w = oper.mul_ApE(eqn, opts, 'N', pc, 'N', V(:, j), 'N') ...
                        + eqn.U * (eqn.V' * V(:, j));
                else
                    w = oper.mul_A(eqn, opts, 'N', V(:, j), 'N') + eqn.U * (eqn.V' * V(:, j));
                end
            else
                if bdf
                    w = (opts.bdf.tau * opts.bdf.beta) * ...
                        oper.mul_ApE(eqn, opts, 'N', pc, 'N', V(:, j), 'N');
                else
                    w = oper.mul_A(eqn, opts, 'N', V(:, j), 'N');
                end
            end
        end
    end
    
    
%     b0 = w;
    for k=1:2 %repeated MGS
        for i = 1 : j
            g = V(:, i)' * w;
            H(i, j) = H(i, j) + g;
            w = w -   V(:, i) * g; 
        end
    end
    beta = norm(w);
    H(j + 1, j) = beta;
    
end

V(:, k + 1) = (1.0 / beta) * w;

end

