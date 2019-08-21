function LQR_FDM_unstable(n0, n_unst, istest)
% Computes a stabilizing feedbacks for primal and dual systems by 
% solving the algebraic Riccati equations with both the ZZ' and LDL'
% approximations of the solutions. 
%
% Inputs:
%
% n0          n0^2 gives the dimension of the original model, i.e. n0 is
%             the number of degrees of freedom, i.e. grid points, per
%             spatial direction 
%             (optional; defaults to 20)
%
% n_unst      number of unstable shifts by construction
%             (optional; defaults to 5)
%
% istest      flag to determine whether this demo runs as a CI test or 
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
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
narginchk(0,3);
if nargin<1, n0 = 20; end
if nargin<2, n_unst = 5; end
if nargin<3, istest = 0; end

%% construct unstable matrix
% first generate the stable part
A0 = fdm_2d_matrix(n0, '10*x', '1000*y', '0'); 
n_stable = size(A0, 1);
n_unstable = n_unst; % no. of unstable eigenvalues
n = n_stable+n_unstable;
B = ones(n, 5);

% now add the instability by constructing a known antistable part and
% corresponding initial stabilizing feedback
Y = eye(n_unstable); % fix "projected feedback" (the Bernoulli/Lyap. sol)

% grab the B part corresponding to the antistable subsystem
Bt = B(n_stable+1 : n, : );
% initial feedback mirrors unstable eigenvals to negative halfplane, these
% may be computed as shifts and make (A+pI) singular  
Vl = [zeros(n_stable,n_unstable); eye(n_unstable)];
K0 = Vl*Y*Bt; 
if exist('lyap','file')% extra antistable part of A
    Au = lyap(Y, -Bt * Bt'); 
else
    Au = lyap2solve(Y, (-Bt * Bt')');
end


%set the full A and corresponding C together with the initial feedback
eqn.A_ = blkdiag(A0,Au);
opts.nm.K0 = K0';
C = ones(10, n);

eqn.haveE = 0;

% set operation
oper = operatormanager('default');

%% global options
opts.norm = 'fro';
%% ADI tolerances and maximum iteration number
opts.adi.maxiter = 200;
opts.adi.res_tol = 1e-10;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;
opts.adi.accumulateDeltaK = 1;
opts.adi.accumulateK = 0;
opts.adi.compute_sol_fac = 1;

%% shift parameters via projection
opts.shifts.num_desired = 5;
opts.shifts.method = 'projection';

%% Newton tolerances and maximum iteration number
opts.nm.maxiter = 15;
opts.nm.res_tol = 1e-8;
opts.nm.rel_diff_tol = 1e-16;
opts.nm.info = 1;
opts.nm.linesearch = 1;
opts.nm.accumulateRes = 1;
opts.nm.res=struct('maxiter',10,'tol',1e-6,'orth',0);

%% for shift banned eigenvalues
opts.shifts.banned = -eig(Au);
opts.shifts.banned_tol = 1e-6;

%% solve Riccati equation
for type=['T','N']
    eqn.type=type;
    eqn.B = B;
    eqn.C = C;
    tic;
    [outB, eqn1, opts1, oper1]=mess_lrnm(eqn, opts, oper);
    toc;
    
    disp('size outB.Z:');
    disp(size(outB.Z));
    
    
    %% check residual norm
    if eqn.type == 'T'
        res1 = mess_res2_norms(outB.Z,'riccati',eqn1,opts1,oper1,opts1.nm,[]) ...
            / norm(eqn.C*eqn.C');
        res2 = abs(eigs(@(x) eqn.A_'*(outB.Z*((outB.Z'*x)))+(outB.Z*((outB.Z'*(eqn.A_*x))))...
            +eqn.C'*(eqn.C*x)-(outB.K)'*((outB.K)*x),n,1,'LM'))...
            /norm(eqn.C*eqn.C');
    else
        res1 = mess_res2_norms(outB.Z,'riccati',eqn1,opts1,oper1,opts1.nm,[]) ...
            / norm(eqn.B'*eqn.B);
        res2 = abs(eigs(@(x) eqn.A_*(outB.Z*((outB.Z'*x)))+(outB.Z*((outB.Z'*(eqn.A_'*x))))...
            +eqn.B*(eqn.B'*x)-(outB.K)'*((outB.K)*x),n,1,'LM'))...
            /norm(eqn.B'*eqn.B);
    end
    fprintf('Newton: %e \t mess_res2_norms: %e \t eigs: %e \n', ...
        [outB.res(end), res1, res2]);
    
    %% print output
    if istest
        if min(outB.res)>=opts.nm.res_tol
            error('MESS:TEST:accuracy','unexpectedly innacurate result');
        end
    else
        figure(1);
        semilogy(outB.res);
        title('0= C^TC + A^TXE + E^TXA -E^TXBB^TXE');
        xlabel('number of iterations');
        ylabel('normalized residual norm');
        pause(1);
    end
    %% Test with LDL_T formulation
    fprintf('\n\n');
    disp('Repeat computations with LDL^T representation of data and solution');
    opts.LDL_T = 1;
    if eqn.type == 'T'
        eqn.S = diag([4,4,9,9,16,16,64,64,81,81]);
        eqn.C = diag([4,4,9,9,16,16,64,64,81,81].^(-0.5)) * eqn.C;
    else
        eqn.S = diag([4,4,9,9,16]);
        eqn.B = eqn.B *diag([4,4,9,9,16].^(-0.5));
    end
    %% solve Riccati equation
    tic;
    [outB2, eqn2, opts2, oper2]=mess_lrnm(eqn, opts, oper);
    toc;
    
    %% print output
    if istest
        if min(outB2.res)>=opts.nm.res_tol
            error('MESS:TEST:accuracy','unexpectedly innacurate result');
        end
    else
        figure(2);
        semilogy(outB2.res);
        title('0= C^TC + A^TXM + M^TXA -M^TXBB^TXM');
        xlabel('number of iterations');
        ylabel('normalized residual norm');
        pause(1);
    end
    disp('size outB2.Z:');
    disp(size(outB2.Z));
    
    
    %% check residual norm
    if eqn.type == 'T'
        res0 = norm(eqn.C'* eqn.S * eqn.C);
        eqn2.S_diag = diag(outB2.S);
        res1 = mess_res2_norms(outB2.Z,'riccati',eqn2,opts2,oper2,opts1.nm,outB2.D)/res0;
        outB2.ZD = outB2.Z*outB2.D;
        res2 = abs(eigs(@(x) eqn.A_'*(outB2.ZD*((outB2.Z'*x)))+(outB2.ZD*((outB2.Z'*(eqn.A_*x))))...
               +eqn.C'*eqn.S*(eqn.C*x)-(outB2.K)'*((outB2.K)*x),n,1,'LM'))...
               /res0;
    else
        res0 = norm((eqn.B * eqn.S) * eqn.B');
        eqn2.S_diag = diag(outB2.S);
        res1 = mess_res2_norms(outB2.Z,'riccati',eqn2,opts2,oper2,opts1.nm,outB2.D)/res0;
        outB2.ZD = outB2.Z*outB2.D;
        res2 = abs(eigs(@(x) eqn.A_*(outB2.ZD*((outB2.Z'*x)))+(outB2.ZD*((outB2.Z'*(eqn.A_'*x))))...
               +eqn.B*eqn.S*(eqn.B'*x)-(outB2.K)'*((outB2.K)*x),n,1,'LM'))...
               /res0;
    end
    fprintf('inner: %e \t mess_res2_norms: %e \t eigs: %e \n', ...
        [outB2.res(end), res1, res2]);
    % only needed to make the next iteration work with ZZ' again
    if opts.LDL_T, opts=rmfield(opts,'LDL_T'); end 
end