function [out, eqn, opts, oper] = mess_lrri(eqn, opts, oper)
%% function [out, eqn, opts, oper] = mess_lrri(eqn, opts, oper)
% Solve continuous-time algebraic H-infinity Riccati equations
%
%   (N) A*X*E' + E*X*A' + E*X*(C1'*C1 - C2'*C2)*X*E' + B1*B1' = 0
%   (T) A'*X*E + E'*X*A + E'*X*(B1*B1' - B2*B2')*X*E + C1'*C1 = 0
%
% for a full-rank factor of X = Z*Z' via the Riccati iteration.
% Using directly and indirectly Matlab M.E.S.S. structures and functions
% for solving definite Riccati equations.
%
% Matrix A can have the form A = Ã + U*V' if U (eqn.U) and V (eqn.V) are
% provided U and V are dense (n x m5) matrices and shoud satisfy m5 << n
%
% Input/Output
%   eqn                 struct contains data for equations
%
%   opts                struct contains parameters for the algorithm
%
%   oper                struct contains function handles for operation
%                       with A and E
%
% Output
%   out                 struct containing solutuions and output information
%
% Input fields in struct eqn:
%   eqn.A_      sparse (n x n) matrix A
%
%   eqn.E_      sparse (n x n) matrix E
%
%   eqn.B1      dense (n x m1) matrix B1
%
%   eqn.B2      dense (n x m2) matrix B2
%               (only if eqn.type = 'T')
%
%   eqn.C1      dense (m3 x n) matrix C1
%
%   eqn.C2      dense (m4 x n) matrix C2
%               (only if eqn.type = 'N')
%
%   eqn.U       dense (n x m5) matrix U
%               (optional, required if eqn.V is present)
%
%   eqn.V       dense (n x m5) matrix V
%               (optional, required if eqn.U is present)
%
%   eqn.type    possible  values: 'N', 'T'
%               determining whether (N) or (T) is solved
%               (optional, default 'N')
%
%   eqn.haveE   possible  values: 0, 1, false, true
%               if haveE = 0: matrix E in eqn.E_ is assumed to be identity
%               (optional, default 0)
%
%   eqn.haveUV  possible  values: 0, 1, false, true
%               if haveUV = 1: U = [U1, U2] and V = [V1, V2]
%               if K or DeltaK are accumulated during the iteration they
%               use only U2 and V2. U1 and V1 can be used for an external
%               rank-k update of the operator.
%               The size of U1 and V1 can be given via eqn.sizeUV1.
%               (optional, default: 0)
%
%   eqn.sizeUV1 possible values: nonnegative integer
%               if a stabilizing feedback is given via U = [U1, U2] and
%               V = [V1, V2] in U2 or V2, eqn.widthU1 indicates how
%               many beginning columns of U and V does not be
%               (optional, default: size(eqn.U, 2))
%
%
% Input fields in struct opts:
%
%   opts.norm               possible  values: 2, 'fro'
%                           use 2-norm (2) or Frobenius norm ('fro') to
%                           compute residual and relative change norms
%                           (optional, default: 'fro')
%
%   opts.ri.Z0              possible values: dense (n x m5) matrix
%                           initial solution factor from the
%                           corresponding LQG problem
%                           (optional, default: [])
%
%   opts.ri.res_tol          possible  values: scalar >= 0
%                           stopping tolerance for the relative
%                           Riccati residual norm; if res_tol = 0 the
%                           relative residual norm is not evaluated
%                           (optional, default: 0)
%
%   opts.ri.rel_diff_tol           possible  values: scalar >= 0
%                           stopping tolerance for the relative
%                           change of the Riccati solution Z;
%                           if res_tol = 0 the relative
%                           change is not evaluated
%                           (optional, default: 0)
%
%   opts.ri.trunc_tol     possible values: scalar >= 0
%                           relative accuracy for column compression
%                           with singular value decomposition
%                           (optional, default: 0)
%
%   opts.ri.riccati_solver  possible values: 'radi', 'newton'
%                           choice of inner Riccati equation solver,
%                           further optional parameters depend on this
%                           'radi'   -> see mess_lrradi
%                           'newton' -> see mess_lrnm
%                           (Note that depending on the initialization,
%                           the first Riccati iteration step might be
%                           performed with the opposite solver.)
%                           (optional, default: 'radi')
%
%   opts.ri.lqg_solver      possible values: 'radi', 'newton'
%                           choice of the Riccati equation solver for
%                           the LQG problem in the first iteration step
%                           'radi'   -> see mess_lrradi
%                           'newton' -> see mess_lrnm
%                           Note: This option is unused if an initial
%                           opts.ri.Z0 is given.
%                           (optional, default: opts.ri.riccati_solver)
%
%   opts.ri.maxiter         possible  values: integer > 0
%                           maximum Riccati iteration number
%                           (optional, default: 10)
%
%   opts.ri.info            possible  values: 0, 1, false, true
%                           turn on (1) or off (0) the status output in
%                           every Riccati iteration step
%                           (optional, default: 0)
%
%   opts.ri.store_lqg       possible values: 0, 1, false, true
%                           if turned on (1) the solution of the LQG
%                           Riccati equation is stored in out.Z_LQG
%                           and the corresponding feedback in out.K_LQG
%                           (optional, default: 0)
%
%   opts.ri.store_solfac    possible values: 0, 1, false, true
%                           if turned on (1) the solution factors
%                           computed by the Riccati equation solvers
%                           are stored in the out.nm and out.radi
%                           structures, otherwise only the information
%                           about the iteration are stored
%                           (optional, default: 0)
%
%   opts.ri.trunc_tol       possible values: scalar > 0
%                           tolerance for rank truncation of the
%                           low-rank solutions (aka column compression)
%                           (optional, default: eps*n)
%
%   opts.ri.trunc_info      possible values: 0, 1, false, true
%                           verbose mode for column compression
%                           (optional, default: 0)
%
%
% If important optional input arguments are missing they may be set to
% default values and a 'MESS:control_data' warning is printed.
% To turn warnings off use warning('OFF', 'MESS:control_data').
%
%
% Output fields in struct out:
%
%   out.Z           low rank solution factor, the solution is X = Z*Z'
%
%   out.K           stabilizing feedback matrix
%
%   out.Z_LQG       low rank solution factor of the corresponding LQG
%                   problem
%                   (opts.ri.store_lqg = 1)
%
%   out.K_LQG       stabilizing feedback matrix of the corresponding LQG
%                   problem
%                   (opts.ri.store_lqg = 1)
%
%   out.niter       number of Riccati iteration steps
%
%   out.res         array of relative Riccati iteration residual norms
%                   (opts.ri.res_tol ~= 0)
%
%   out.rc          array of relative Riccati iteration change norms
%                   (opts.ri.rel_diff_tol ~= 0)
%
%   out.res0        norm of the normalization residual term
%
%   out.nm          struct with output of all Newton iterations
%
%   out.radi        struct with output of all RADI iterations
%
%
%   uses operatorfunctions init and mul_E, mul_E_pre, mul_E_post
%   and further indirectly in the inner Riccati solver
%
%   See also mess_lrnm, mess_lradi, mess_para,
%   mess_galerkin_projection_acceleration, operatormanager.

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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check system data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'haveE')) || isempty(eqn.haveE)
    eqn.haveE = 0;
end

% Initialize function operator.
[result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A', 'E');
if not(result)
    error( 'MESS:control_data', ...
        'system data is not completely defined or corrupted' );
end

% Check type of equation.
if not(isfield(eqn, 'type'))
    eqn.type = 'N';
    warning('MESS:control_data',['Unable to determine type of equation.'...
        'Falling back to type ''N''']);
elseif (eqn.type ~= 'N') && (eqn.type ~= 'T')
    error('MESS:equation_type', ...
        'Equation type must be either ''T'' or ''N''');
end

% make sure the corresponding matrices from quadratic term are well
% defined and the first right hand side is dense so that the resulting
% factor is densly stored.
if not(isfield(eqn, 'B1')) || not(isnumeric(eqn.B1))
    error('MESS:control_data', 'eqn.B1 is not defined or corrupted');
end

if not(isfield(eqn, 'C1')) || not(isnumeric(eqn.C1))
    error( 'MESS:control_data', 'eqn.C1 is not defined or corrupted');
end

if issparse(eqn.B1), eqn.B1 = full(eqn.B1); end
if issparse(eqn.C1), eqn.C1 = full(eqn.C1); end

if eqn.type == 'T'
    if not(isfield(eqn, 'B2')) || not(isnumeric(eqn.B2))
        error('MESS:control_data', 'eqn.B2 is not defined or corrupted');
    end
    
    if issparse(eqn.B2), eqn.B2 = full(eqn.B2); end
else
    if not(isfield(eqn, 'C2')) || not(isnumeric(eqn.C2))
        error('MESS:control_data', 'eqn.C2 is not defined or corrupted');
    end
    
    if issparse(eqn.C2), eqn.C2 = full(eqn.C2); end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rank-k update system data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(eqn, 'haveUV')) || isempty(eqn.haveUV) || not(eqn.haveUV)
    eqn.haveUV  = 0;
    eqn.sizeUV1 = 0;
    eqn.U       = [];
    eqn.V       = [];
else
    if isnumeric(eqn.U) ...
            && isnumeric(eqn.V) ...
            && size(eqn.U, 1) == size(eqn.V, 1) ...
            && size(eqn.U, 2) == size(eqn.V, 2)
        
        if issparse(eqn.V), eqn.V = full(eqn.V); end
        if issparse(eqn.U), eqn.U = full(eqn.U); end
    else
        error('MESS:control_data', ...
            'Inappropriate data of low rank updated operator (eqn.U and eqn.V)');
    end
end

% Check for size of constant term in U and V.
if eqn.haveUV
    if not(isfield(eqn, 'sizeUV1')) || isempty(eqn.sizeUV1)
        eqn.sizeUV1 = size(eqn.U, 2);
    else
        assert(isnumeric(eqn.sizeUV1) ...
            && (eqn.sizeUV1 <= size(eqn.U, 2)), ...
            'MESS:control_data', ...
            'Inappropriate size of low rank updated operator (eqn.U and eqn.V)');
    end
end
init_sizeUV1 = eqn.sizeUV1;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize required usf for multiplications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eqn.haveE
    [eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for Riccati Iteration control structure in options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts, 'ri')) || not(isstruct(opts.ri))
    error('MESS:control_data', ['No ri control data found in ', ...
        'options structure.']);
end

if not(isfield(opts.ri, 'maxiter')) || not(isnumeric(opts.ri.maxiter))
    warning( 'MESS:control_data', ...
        ['Missing or corrupted maxiter field. ', ...
        'Switching to default opts.ri.maxiter = 10.']);
    opts.ri.maxiter = 10;
end

if not(isfield(opts.ri, 'res_tol')) || not(isnumeric(opts.ri.res_tol))
    warning('MESS:control_data', ...
        ['Missing or corrupted res_tol field. ', ...
        'Switching to default opts.ri.res_tol = 0.']);
    opts.ri.res_tol = 0;
end

if not(isfield(opts.ri, 'rel_diff_tol')) || not(isnumeric(opts.ri.rel_diff_tol))
    warning('MESS:control_data', ...
        ['Missing or corrupted rel_diff_tol field. ', ...
        'Switching to default opts.ri.rel_diff_tol = 0.']);
    opts.ri.rel_diff_tol = 0;
end

if not(isfield(opts.ri, 'compres_tol')) || not(isnumeric(opts.ri.compres_tol))
    warning('MESS:control_data', ...
        ['Missing or corrupted compres_tol field. ', ...
        'Switching to default opts.ri.compres_tol = 0.']);
    opts.ri.compres_tol = 0;
end

if not(isfield(opts.ri, 'riccati_solver')) || isempty(opts.ri.riccati_solver)
    warning('MESS:control_data', ...
        ['Missing or corrupted riccati_solver field. ', ...
        'Switching to default opts.ri.riccati_solver = ''radi''.']);
    
    opts.ri.riccati_solver    = 'radi';
    riccati_solver            = @mess_lrradi;
    opts.radi.compute_sol_fac = 1;
elseif strcmpi(opts.ri.riccati_solver, 'radi')
    riccati_solver            = @mess_lrradi;
    opts.radi.compute_sol_fac = 1;
elseif strcmpi(opts.ri.riccati_solver, 'newton')
    riccati_solver           = @mess_lrnm;
    opts.adi.compute_sol_fac = 1;
else
    error('MESS:notimplemented', ...
        'The requested Riccati solver is not implemented.');
end

if not(isfield(opts.ri, 'lqg_solver')) || isempty(opts.ri.lqg_solver)
    lqg_solver = riccati_solver;
elseif strcmpi(opts.ri.lqg_solver, 'radi')
    lqg_solver                = @mess_lrradi;
    opts.radi.compute_sol_fac = 1;
elseif strcmpi(opts.ri.lqg_solver, 'newton')
    lqg_solver               = @mess_lrnm;
    opts.adi.compute_sol_fac = 1;
else
    error('MESS:notimplemented', ...
        'The requested Riccati solver is not implemented.');
end

if not(isfield(opts.ri, 'info'))
    opts.ri.info = 0;
else
    if not(isnumeric(opts.ri.info)) && not(islogical(opts.ri.info))
        error( 'MESS:control_data', ...
            'opts.ri.info parameter must be logical or numeric.');
    end
end

if not(isfield(opts.ri, 'store_lqg')) || isempty(opts.ri.store_lqg)
    opts.ri.store_lqg = 0;
else
    if not(isnumeric(opts.ri.store_lqg)) && not(islogical(opts.ri.store_lqg))
        error( 'MESS:control_data', ...
            'opts.ri.store_lqg parameter must be logical or numeric.');
    end
end

if not(isfield(opts.ri, 'store_solfac')) || isempty(opts.ri.store_solfac)
    opts.ri.store_solfac = 0;
else
    if not(isnumeric(opts.ri.store_solfac)) && not(islogical(opts.ri.store_solfac))
        error( 'MESS:control_data', ...
            'opts.ri.store_solfac parameter must be logical or numeric.');
    end
end

% Check for residual norm.
if not(isfield(opts, 'norm')) ...
        || (not(strcmp(opts.norm, 'fro')) && ...
        (not(isnumeric(opts.norm)) || opts.norm ~= 2))
    warning('MESS:control_data', ...
        ['Missing or Corrupted opts.norm field.', ...
        'Switching to default: ''fro''']);
    opts.norm = 'fro';
end

% Check for incompatible shift selection.
ham_shifts = 0;
if strcmpi(func2str(riccati_solver), 'radi') ...
        && strcmpi(func2str(lqg_solver), 'newton') ...
        && isfield(opts, 'shifts') ...
        && isfield(opts.shifts, 'method') ...
        && strcmpi(opts.shifts.method, 'gen-ham-opti')
    warning('MESS:control_data', ...
        ['The chosen shift method is not usable in the LQG step. ', ...
         'The shift method will be changed for this step to ', ...
         '''projection'' and for the inner iteration back to its ', ...
         'original state.']);
    
    ham_shifts         = 1;
    opts.shifts.method = 'projection';
end

if not(isfield(opts.ri,'trunc_tol')),  ...
    opts.ri.trunc_tol=eps * oper.size(eqn, opts); end
if not(isfield(opts.ri, 'trunc_info')), opts.ri.trunc_info = 0; end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List all currently unsupported options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(opts, 'LDL_T') && opts.LDL_T
    error('MESS:notimplemented', ...
        'The LDL_T factorization type is not supported in this function.');
end
opts.LDL_T = false; % We need this to apply oper.init_res later.

if isfield(opts, 'bdf') && not(isempty(opts.bdf))
    error( 'MESS:control_data', 'Options bdf not supported.');
end

if isfield(opts, 'rosenbrock') && not(isempty(opts.rosenbrock))
    error( 'MESS:control_data', 'Options rosenbrock not supported.');
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All checks done. Here comes the real work!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = [];

if opts.ri.res_tol
    res = zeros(1, opts.ri.maxiter);
else
    res = [];
end

if opts.ri.rel_diff_tol
    rc = zeros(1, opts.ri.maxiter);
else
    rc = [];
end

if opts.ri.rel_diff_tol
    normZ = 0;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATION PHASE 1: Solve the LQG problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.ri, 'Z0')) || isempty(opts.ri.Z0)
    % Setup MESS equation structure for initial Riccati equation.
    if eqn.type == 'T'
        [~, res0, eqn, opts, oper] = ...
            oper.init_res(eqn, opts, oper, eqn.C1');
        eqn.B = eqn.B2;
        eqn.C = eqn.C1;
    else
        [~, res0, eqn, opts, oper] = ...
            oper.init_res(eqn, opts, oper, eqn.B1);
        eqn.C = eqn.C2;
        eqn.B = eqn.B1;
    end
    
    [out_riccati, eqn, opts, oper] = lqg_solver(eqn, opts, oper);
    
    % Store the solution of the LQG problem.
    if opts.ri.store_lqg
        out.Z_LQG = out_riccati.Z;
        out.K_LQG = out_riccati.K;
    end
    
    % Store information about Riccati equation solver.
    tmp = out_riccati;
    if not(opts.ri.store_solfac)
        if isfield(tmp, 'Z'), tmp = rmfield(tmp, 'Z'); end
        if isfield(tmp, 'K'), tmp = rmfield(tmp, 'K'); end
        if isfield(tmp, 'res_fact'), tmp = rmfield(tmp, 'res_fact'); end
    end
    switch func2str(lqg_solver)
        case 'mess_lrradi'
            out.radi(1) = tmp;
        case 'mess_lrnm'
            out.nm(1) = tmp;
    end

    if ham_shifts % *** REMOVE IN LATER VERSION: MAKE SHIFTS COMPATIBLE
        opts.shifts.method = 'gen-ham-opti';
    end
else
    % Compute zero solution residual for normalization.
    if eqn.type == 'T'
        [~, res0, eqn, opts, oper] = ...
            oper.init_res(eqn, opts, oper, eqn.C1');
        eqn.B = eqn.B2;
    else
        [~, res0, eqn, opts, oper] = ...
            oper.init_res(eqn, opts, oper, eqn.B1);
        eqn.C = eqn.C2;
    end
    
    % Take the given initial solution.
    out_riccati.Z = opts.ri.Z0;
end

% Remove initial feedback for Newton method.
isnmK0 = 0;
if isfield(opts, 'nm') && isfield(opts.nm, 'K0')
    isnmK0  = 1;
    nmK0    = opts.nm.K0;
    opts.nm = rmfield(opts.nm, 'K0');
end

% Remove initial matrices for RADI method.
isradiZ0 = 0;
isradiY0 = 0;
isradiK0 = 0;
isradiW0 = 0;
if isfield(opts, 'radi')
    if isfield(opts, 'radi') && isfield(opts.radi, 'Z0')
        isradiZ0  = 1;
        radiZ0    = opts.radi.Z0;
        opts.radi = rmfield(opts.radi, 'Z0');
    end
    
    if isfield(opts.radi, 'Y0')
        isradiY0  = 1;
        radiY0    = opts.radi.Y0;
        opts.radi = rmfield(opts.radi, 'Y0');
    end
    
    if isfield(opts.radi, 'K0')
        isradiK0  = 1;
        radiK0    = opts.radi.K0;
        opts.radi = rmfield(opts.radi, 'K0');
    end
    
    if isfield(opts.radi, 'W0')
        isradiW0  = 1;
        radiW0    = opts.radi.W0;
        opts.radi = rmfield(opts.radi, 'W0');
    end
end

% Check if the shift history needs to be reset for RADI.
if isfield(opts, 'shifts') ...
        && isfield(opts.shifts, 'method') ...
        && strcmpi(opts.shifts.method, 'gen-ham-opti') ...
        && isfield(opts.shifts, 'history') ...
        && exist('W', 'var')
    if eqn.type == 'T'
        if mod(opts.shifts.history, size(eqn.B1, 2))
            k = ceil(opts.shifts.history / size(W, 2));
            opts.shifts.history = k * size(eqn.B1, 2);
            warning('MESS:control_data', ...
                ['Size of the residual changed after LQG problem. ', ...
                'The parameter opts.shifts.history is reset to %d.'], ...
                opts.shifts.history);
        end
    else
        if mod(opts.shifts.history, size(eqn.C1, 1))
            k = ceil(opts.shifts.history / size(W, 2));
            opts.shifts.history = k * size(eqn.C1, 1);
            warning('MESS:control_data', ...
                ['Size of the residual changed after LQG problem. ', ...
                'The parameter opts.shifts.history is reset to %d.'], ...
                opts.shifts.history);
        end
    end
end

% Initialize storage for the computed feedback.
if eqn.type == 'T'
    [n, m1] = size(eqn.B1);
    m2      = size(eqn.B2, 2);
    m12     = m1 + m2;
    eqn.U   = [eqn.U(:, 1:eqn.sizeUV1), eqn.B1, -eqn.B2];
    eqn.V   = [eqn.V(:, 1:eqn.sizeUV1), zeros(n, m12)];
else
    [m1, n] = size(eqn.C1);
    m2      = size(eqn.C2, 1);
    m12     = m1 + m2;
    eqn.U   = [eqn.U(:, 1:eqn.sizeUV1), zeros(n, m12)];
    eqn.V   = [eqn.V(:, 1:eqn.sizeUV1), eqn.C1', -eqn.C2'];
end
eqn.haveUV  = 1;
eqn.sizeUV1 = size(eqn.V, 2);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERSTION PHASE 2: Riccati Iteration (solve the residual equations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:1:opts.ri.maxiter
    % Accumulate new solution.
    if (opts.ri.trunc_tol > 0) && (j > 1)
        Z = mess_column_compression([Z, out_riccati.Z], 'N', [], ...
            opts.ri.trunc_tol, opts.ri.trunc_info);

    else
        Z(:, end+1:end+size(out_riccati.Z, 2)) = out_riccati.Z;
    end
    
    % Update the constant term and error variables.
    if eqn.type == 'T'
        if eqn.haveE
            eqn.C = (oper.mul_E( eqn, opts, 'T', out_riccati.Z, 'N' ) ...
                * (out_riccati.Z' * eqn.B1))';
        else
            eqn.C = (eqn.B1' * out_riccati.Z) * out_riccati.Z';
        end
    else
        if eqn.haveE
            eqn.B = oper.mul_E(eqn, opts, 'N', out_riccati.Z, 'N') ...
                * (eqn.C1 * out_riccati.Z)';
        else
            eqn.B = out_riccati.Z * (eqn.C1 * out_riccati.Z)';
        end
    end
    
    % Set the next rank-k update
    if eqn.type == 'T'
        if eqn.haveE
            eqn.V(:, end-m12+1:end) = ...
                oper.mul_E(eqn, opts, 'T', Z, 'N') ...
                * (Z' * [eqn.B1, eqn.B2]);
        else
            eqn.V(:, end-m12+1:end) = Z * (Z' * [eqn.B1, eqn.B2]);
        end
    else
        if eqn.haveE
            eqn.U(:, end-m12+1:end) = ...
                oper.mul_E(eqn, opts, 'N', Z, 'N') ...
                * (Z' * [eqn.C1', eqn.C2']);
        else
            eqn.U(:, end-m12+1:end) = Z * (Z' * [eqn.C1', eqn.C2']);
        end
    end
    
    % Compute convergence measures.
    if opts.ri.res_tol
        if eqn.type == 'T'
            res(j) = norm(eqn.C * eqn.C', opts.norm) / res0;
        else
            res(j) = norm(eqn.B' * eqn.B, opts.norm) / res0;
        end
    end
    
    if opts.ri.rel_diff_tol
        normY = sum(sum(out_riccati.Z.^2));
        normZ = normZ + normY;
        rc(j) = normY / normZ;
    end
    
    % Print status information.
    if opts.ri.info
        if opts.ri.rel_diff_tol && opts.ri.res_tol
            fprintf(1, ...
                ['RI step: %4d  normalized residual: %e ' ...
                'relative change in Z: %e\n'], ...
                j, res(j), rc(j));
        elseif opts.ri.res_tol
            fprintf(1, 'RI step: %4d normalized residual: %e\n', ...
                j, res(j));
        elseif opts.ri.rel_diff_tol
            fprintf(1, 'RI step: %4d relative change in Z: %e\n', ...
                j, rc(j));
        end
        
        if isfield(out_riccati, 'adi')
            fprintf(1, '               number of Newton steps: %4d\n\n', ...
                out_riccati.niter);
        elseif isfield(out_riccati, 'niter')
            fprintf(1, '               number of RADI steps: %4d\n\n', ...
                out_riccati.niter);
        end
    end
    
    % Evaluate stopping criteria.
    if (opts.ri.res_tol && (res(j) < opts.ri.res_tol)) ...
            || (opts.ri.rel_diff_tol && (rc(j) < opts.ri.rel_diff_tol)) ...
            || (j >= opts.ri.maxiter) 
        break;
    end
    
    % Solve the next residual equation.
    [out_riccati, eqn, opts, oper] = riccati_solver(eqn, opts, oper);
    
    % Store information about Riccati equation solver.
    tmp = out_riccati;
    if not(opts.ri.store_solfac)
        if isfield(tmp, 'Z'), tmp = rmfield(tmp, 'Z'); end
        if isfield(tmp, 'K'), tmp = rmfield(tmp, 'K'); end
        if isfield(tmp, 'res_fact'), tmp = rmfield(tmp, 'res_fact'); end
    end
    switch func2str(riccati_solver)
        case 'mess_lrradi'
            out.radi(j+1) = tmp;
        case 'mess_lrnm'
            out.nm(j+1) = tmp;
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.Z = Z;

out.niter = j;

if eqn.type == 'T'
    out.K = eqn.V(:, end-m12+1:end)';
else
    out.K = eqn.U(:, end-m12+1:end)';
end

if opts.ri.res_tol
    out.res = res(1:j);
end

if opts.ri.rel_diff_tol
    out.rc = rc(1:j);
end

out.res0 = res0;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eqn.sizeUV1 = init_sizeUV1;
if (size(eqn.V, 2) > eqn.sizeUV1) || (size(eqn.U, 2) > eqn.sizeUV1)
    % Cut off the stabilizing feedback.
    eqn.V = eqn.V(:, 1:eqn.sizeUV1);
    eqn.U = eqn.U(:, 1:eqn.sizeUV1);
end

if isempty(eqn.V) || isempty(eqn.U)
    % Enforce empty matrices and parameters.
    eqn.U       = [];
    eqn.V       = [];
    eqn.haveUV  = 0;
    eqn.sizeUV1 = 0;
end

% Delete overwritten right hand-side and quadratic term.
eqn = rmfield(eqn, 'B');
eqn = rmfield(eqn, 'C');

% Rebuild initial values in option struct.
if isnmK0, opts.nm.K0 = nmK0; end
if isradiZ0, opts.radi.Z0 = radiZ0; end
if isradiY0, opts.radi.Y0 = radiY0; end
if isradiK0, opts.radi.K0 = radiK0; end
if isradiW0, opts.radi.W0 = radiW0; end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalize required usfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eqn.haveE
    [eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
end
