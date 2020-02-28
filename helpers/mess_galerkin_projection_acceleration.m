function [Z, D, S]=mess_galerkin_projection_acceleration(Z, type, ...
                                                  eqn, oper, fopts, D)
%  Galerkin projection onto subspace spanned by low-rank factor Z of ADI
%  method for solving FXE^T + EXF^T = -GG^T.
%
% Input:
%  Z                        Low-rank factor Z
%
%  type                     possible values: 'LE','CARE'
%                           determines whether a Lyapunov ('LE') or a
%                           Riccati ('CARE') equation should be
%                           projected
%
%  eqn                      structure with data for A, E, G
%                           eqn.E(optional, eqn.haveE specifies whether it
%                           is there) in the above equation with ZZ'
%                           approximating X
%
%  oper                     structure contains function handles for
%                           operations with A, E
%
%  fopts                    options structure that should contain following
%                           members
%
%  fopts.adi                options structure for ADI method
%
%  fopts.nm                 options structure for Newton method
%
%  xopts                    fopts.adi or fopts.nm depending on type
%
%  xopts.projection.ortho   possible  values: 0, 1, false, true
%                           implicit (0) or explicit (1)
%                           orthogonalization of Z in NM
%                           i.e., explicit orthogonalization via orth().
%                           (optional, default: 1)
%
%  xopts.projection.meth    method for solving projected Lyapunov or
%                           Riccati equation xopts is fopts.adi or fopts.nm
%                           depending on type possible  values: 'lyapchol',
%                           'lyap_sgn_fac','lyap','lyapunov','lyap2solve', 
%                           'care','care_nwt_fac','mess_dense_nm'
%                           (optional, default: 'lyap' or 'care' depending
%                           on type)
%
%  D                        solution factor D for LDL^T formulation
%
% Output:
%  Z                        Updated solution factor Z after prolongation
%
%  D                        Updated solution factor D after prolongation
%
%  S                        Updated solution factor S after prolongation
%
% uses oparatorfunctions mul_A, mul_E, mul_ApE

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

factorize=1;

switch type
    case {'LE'}
        opts=fopts.adi;
    case {'CARE'}
        opts=fopts.nm;
    otherwise
        error('MESS:GP_type',...
            ['type has to be ''LE'' or ''CARE'' in Galerkin projection '...
            'acceleration.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.projection,'ortho')) || isempty(opts.projection.ortho)
    opts.projection.ortho=1;
end
% if opts.projection.meth is not set use default, availability is checked
% later
if not(isfield(opts.projection,'meth')) || isempty(opts.projection.meth)
    set_default = 1;
    switch type
        case 'LE'
            opts.projection.meth='lyap';
        case 'CARE'
            opts.projection.meth='care';
    end
else
    set_default = 0;
end

if strcmp(opts.projection.meth, 'riccati')
    opts.projection.meth = 'care';
end
if strcmp(opts.projection.meth, 'lyapunov')
    opts.projection.meth = 'lyap';
end

if not(exist(opts.projection.meth,'file'))
    switch type
        case 'LE'
            meth_non_ex = opts.projection.meth;
            if exist('lyap', 'file')
                opts.projection.meth='lyap';
            elseif exist('lyapchol', 'file')
                opts.projection.meth='lyapchol';
            elseif exist('lyap_sgn_fac', 'file')
                opts.projection.meth='lyap_sgn_fac';
            elseif exist('lyap2solve', 'file')
                opts.projection.meth='lyap2solve';
            else
                error('MESS:GP_missing_solver',...
                    ['Galerkin projection acceleration was unable to '...
                    'find Lyapunov solver']);
            end
        case 'CARE'
            meth_non_ex = opts.projection.meth;
            if exist('icare', 'file')
                opts.projection.meth='icare';
            elseif exist('care', 'file')
                opts.projection.meth='care';
            elseif exist('care_nwt_fac', 'file')
                opts.projection.meth='care_nwt_fac';
            elseif exist('mess_dense_nm', 'file')
                opts.projection.meth='mess_dense_nm';
            else
                error('MESS:GP_missing_solver',...
                    ['Galerkin projection acceleration was unable to '...
                    'find Riccati solver']);
            end
    end
    if not(set_default)
        warning('MESS:GP_missing_solver',...
            ['Galerkin projection acceleration was unable to find'...
            ' solver ''%s'', switched to ''%s'''], meth_non_ex, ...
            opts.projection.meth);
    end
end

if not(isfield(fopts,'rosenbrock')), fopts.rosenbrock=[]; end
if isstruct(fopts.rosenbrock)&&isfield(fopts.rosenbrock,'tau')
    rosenbrock = 1;
    if fopts.rosenbrock.stage == 1
        pc = -1 / (2 * fopts.rosenbrock.tau);
    else % stage 2
        pc = -1 / (2 * fopts.rosenbrock.tau * fopts.rosenbrock.gamma);
    end
else
    rosenbrock = 0;
end
if not(isfield(fopts,'bdf')), fopts.bdf=[]; end
if isstruct(fopts.bdf) && isfield(fopts.bdf, 'tau') ...
        && isfield(fopts.bdf, 'beta')
    bdf = 1;
    pc = -1 / (2 * fopts.bdf.tau * fopts.bdf.beta);
else
    bdf = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute projector matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fopts.LDL_T
    [Z, D] = mess_column_compression(Z, 'N', D, eps, 0);
    if strcmp(opts.projection.meth, 'care_nwt_fac') ...
            || strcmp(type, 'LE')
        if not(set_default)
            meth_not_supp = opts.projection.meth;
        end
        switch type
            case 'LE'
                if isfield(eqn, 'diagonalized_RHS') && eqn.diagonalized_RHS
                    U = eqn.U_diag;
                else
                    U = 1;
                end
                S = diag(eqn.S_diag);
                if any(diag(S)<0) && ...
                        (strcmp(opts.projection.meth, 'lyapchol') ...
                         || strcmp(opts.projection.meth, 'lyap_sgn_fac'))
                    if exist('lyap', 'file')
                        opts.projection.meth='lyap';
                    elseif exist('lyap2solve', 'file')
                        opts.projection.meth='lyap2solve';
                    else
                        error('MESS:GP_missing_solver',...
                            ['Galerkin projection acceleration was unable to '...
                            'find Lyapunov solver']);
                    end
                end
            case 'CARE'
                [U,S] = eig(eqn.S);
                if any(diag(S)<0)
                    if exist('care', 'file')
                        opts.projection.meth='care';
                    elseif exist('mess_dense_nm', 'file')
                        opts.projection.meth='mess_dense_nm';
                    else
                        error('MESS:GP_missing_solver',...
                            ['Galerkin projection acceleration was unable to '...
                            'find Riccati solver']);
                    end
                end
        end
        if not(set_default) && not(strcmp(meth_not_supp, opts.projection.meth))
            warning('MESS:galerkin_projection_acceleration', ...
                ['Projection method %s not supported for ', ...
                'LDL_T formulation with D indefinite. Switching to %s.'], ...
                meth_not_supp, opts.projection.meth);
        end
    end
elseif opts.projection.ortho
    Z=orth(Z);
else
    [U,S,~]=svd(full(Z'*Z));
    s=diag(S);
    sk=find(s>eps*s(1), 1, 'last' );
    Z=Z*U(:,1:sk)*diag(1./sqrt(s(1:sk)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the projected Matrix equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type
    case 'LE'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The Lyapunov equation case
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lyapunov = 1;
        B=Z'*eqn.G;
        if bdf || rosenbrock
            A = oper.mul_ApE(eqn, fopts,eqn.type,pc,eqn.type,Z,'N');
            if bdf
                A = (fopts.bdf.tau * fopts.bdf.beta) * A;
                if eqn.haveUV
                    if eqn.type=='T'
                        A = A + eqn.V * (eqn.U' * Z);
                    else
                        A = A + eqn.U * (eqn.V' * Z);
                    end
                end
            else % rosenbrock
                if fopts.rosenbrock.stage == 2
                    A = (fopts.rosenbrock.tau * fopts.rosenbrock.gamma) * A;
                end
                if eqn.haveUV
                    if eqn.type=='T'
                        A = A + eqn.V * (eqn.U' * Z);
                    else
                        A = A + eqn.U * (eqn.V' * Z);
                    end
                end
            end
        else
            A = oper.mul_A( eqn, fopts, eqn.type, Z, 'N' );
            if eqn.haveUV
                if eqn.type=='T'
                    A = A + eqn.V * (eqn.U' * Z);
                else
                    A = A + eqn.U * (eqn.V' * Z);
                end
            end
        end
        A = Z' * A;
        if eqn.haveE
            E = Z'*(oper.mul_E(eqn, fopts,eqn.type,Z,'N'));
        else
            E = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Choose solver for the small equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch opts.projection.meth
            case 'lyapchol'
                % in LDL_T with S neg. EV sqrt(S_diag) will be complex
                if fopts.LDL_T
                    B = B * U * sqrt(S);
                    if eqn.haveE
                        XC = lyapchol(A,B,E);
                    else
                        XC = lyapchol(A,B);
                    end
                    [~,S,XC] = svd(XC,'econ');
                    XC = XC';
                    D = S.^2;
                    S = 1;
                else
                    if eqn.haveE
                        XC=lyapchol(A,B,E);
                    else
                        XC=lyapchol(A,B);
                    end
                end
                factorize=0;
                
            case 'lyap_sgn_fac'
                if fopts.LDL_T
                    B = B * U * sqrt(S);
                    XC = lyap_sgn_fac(A',B',E');
                    [~,S,XC] = svd(XC,'econ');
                    XC = XC';
                    D = S.^2;
                    S = 1;
                else
                    XC = lyap_sgn_fac(A',B',E');
                end
                factorize=0;
                
            case {'lyap','lyapunov'}
                if fopts.LDL_T
                    B = B*U*S*U'*B';
                    B = (B + B') / 2; % make sure it's symmetric for lyap
                    if eqn.haveE
                        X = lyap(A,B,[],E);
                    else
                        X = lyap(A,B);
                    end
                else
                    if eqn.haveE
                        X = lyap(A,B*B',[],E);
                    else
                        X = lyap(A,B*B');
                    end
                end
                
            case 'lyap2solve'
                if eqn.haveE
                    EB = E\B;
                    if fopts.LDL_T
                        X = lyap2solve(E\A,EB*U*S*U'*EB');
                    else
                        X = lyap2solve(E\A,EB*EB');
                    end
                else
                    if fopts.LDL_T
                        X = lyap2solve(A,B*U*S*U'*B');
                    else
                        X = lyap2solve(A,B*B');
                    end
                end
        end
    case 'CARE'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The Riccati equation case
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lyapunov = 0;
        if eqn.type=='T'
            opAE = 'N';
            if bdf
                B=Z'*eqn.B * sqrt(fopts.bdf.tau * fopts.bdf.beta);
            else
                B=Z'*eqn.B;
            end
            C=eqn.C*Z;
        else
            opAE = 'T';
            if bdf
                B=Z'*eqn.C' * sqrt(fopts.bdf.tau * fopts.bdf.beta);
            else
                B=Z'*eqn.C';
            end
            C=eqn.B'*Z;
        end
        if bdf
            A = (fopts.bdf.tau * fopts.bdf.beta) * (Z' ...
                * oper.mul_ApE(eqn, fopts,opAE,pc,opAE,Z,'N'));
        else
            A = oper.mul_A( eqn, fopts, opAE, Z, 'N' );
            if eqn.haveUV && eqn.sizeUV1
                if eqn.type=='T'
                    A = A + eqn.U(:, 1:eqn.sizeUV1) ...
                        * (eqn.V(:, 1:eqn.sizeUV1)' * Z);
                else
                    A = A + eqn.V(:, 1:eqn.sizeUV1) ...
                        * (eqn.U(:, 1:eqn.sizeUV1)' * Z);
                end
            end
            A = Z' * A;
        end
        if eqn.haveE
            E = Z'*(oper.mul_E(eqn, fopts,opAE,Z,'N'));
        else
            E = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Choose solver for the small equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch opts.projection.meth
            case {'care', 'riccati'}
                if exist('icare','file')
                    warning('MESS:care',['It seems that icare() is ' ...
                                        'available on your system. ' ...
                                        'We recommend using icare() ' ...
                                        'over using care() following the ' ...
                                        'recommendation by TMW.'])
                end
            
                if fopts.LDL_T
                    X=care(A,B,C'*eqn.S*C,eye(size(B,2)),[],E);
                else
                    X=care(A,B,C'*C,eye(size(B,2)),[],E);
                end
         
            case {'icare'}
                if fopts.LDL_T
                    X=icare(A,B,C'*eqn.S*C,eye(size(B,2)),[],E);
                else
                    X=icare(A,B,C'*C,eye(size(B,2)),[],E);
                end
         
            case 'care_nwt_fac'
                if fopts.LDL_T
                    C = sqrt(S) * U' * C;
                end
                if not(isempty(E))
                    XC = care_nwt_fac([],A/E,B,C/E,1e-12,50);
                else
                    XC = care_nwt_fac([],A,B,C,1e-12,50);
                end
                if fopts.LDL_T
                    [~,S,XC] = svd(XC,'econ');
                    XC = XC';
                    D = S.^2;
                    S = 1;
                end
                factorize=0;
                
            case 'mess_dense_nm'
                    if fopts.LDL_T
                        X=mess_dense_nm(A,B,C,E, [], eqn.S);
                    else
                        X=mess_dense_nm(A,B,C,E);
                    end
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the projected solution was not already computed in factored form
% compute a symmetric factorization now and update the large factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if factorize
    if fopts.LDL_T
        [V, D] = eig(X);
        Z = Z * V;
        S = 1;
    elseif (exist('cholp','file'))
        [XC,P,I]=cholp(X);
        XC=XC*P';
        if I && lyapunov
            warning('MESS:proj_sol_semidef',...
                'The solution of the projected equation was semidefinite.');
        end
    else
        [~,S,V]=svd(X);
        s=diag(S);
        r=find(s>s(1)*eps);
        XC=diag(sqrt(s(r)))*V(:,r)';
    end
end
if exist('XC','var'), Z=Z*XC'; end
