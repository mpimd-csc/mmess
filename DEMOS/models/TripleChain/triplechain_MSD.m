function [M,E,K]=triplechain_MSD(n1,alpha,beta,v)
% function [M,E,K]=triplechain_MSD(n1,alpha,beta,v)
%
% Generates mass spring damper system of three coupled mass sprong damper
% chains as in [1,example 2] with proportional damping. The resulting
% system has dimension 3*n1+1
%
% Output:
%
% M,E,K    mass , damping and stiffness matrices of the system
%
% Input:
%
% n1          length of each of the three coupled chains.
% alpha,beta  coefficients in the proportional damping rule:
%             Dp=alpha*M+beta*K;
% v           viskosity of the three additional dampers
%
% [1] N.Truhar and K.Veselic
%     An efficient method for estimating the optimal dampers' viscosity for
%     linear vibrating systems using Lyapunov equations
%     SIAM J. Matrix Anal. Appl. vol.31 no.1 (2009) pp 18-39
%     https://doi.org/10.1137/070683052.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%



m1=1;
m2=2;
m3=3;
m0=10;

k1=10;
k2=20;
k3=1;
k0=50;

if nargin<2
  alpha=.002;
  beta=alpha;
  v=5e0;
end

M=spdiags([m1*ones(1,n1) m2*ones(1,n1) m3*ones(1,n1) m0]',0,3*n1+1,3*n1+1);
e=ones(n1,1);
K1=spdiags([-e 2*e -e],-1:1,n1,n1);


K=sparse(3*n1+1,3*n1+1);
K(1:n1,1:n1)=k1*K1;
K(n1+1:2*n1,n1+1:2*n1)=k2*K1;
K(2*n1+1:3*n1,2*n1+1:3*n1)=k3*K1;

K(1:n1,end)=-[sparse(1,n1-1) k1]';
K(n1+1:2*n1,end)=-[sparse(1,n1-1) k2]';
K(2*n1+1:3*n1,end)=-[sparse(1,n1-1) k3]';

K(end,1:n1)=-[sparse(1,n1-1) k1];
K(end,n1+1:2*n1)=-[sparse(1,n1-1) k2];
K(end,2*n1+1:3*n1)=-[sparse(1,n1-1) k3];
K(3*n1+1, 3*n1+1)=k1+k2+k3+k0;

%fractional damping is used in [1] but needs full matrices due to the sqrtm
%SM=spdiags([sqrt(m1)*ones(1,n1) sqrt(m2)*ones(1,n1) sqrt(m3)*ones(1,n1)...
%  sqrt(m0)]',0,3*n1+1,3*n1+1);
%E=.02*(2*SM*sqrtm((SM\K)/SM)*SM);

% here we want internal damping based on sparse matrices via Rayleigh damping:
E=alpha*M+beta*K;
E(1,1)=E(1,1)+v;
E(n1,n1)=E(n1,n1)+v;
E(2*n1+1,2*n1+1)=E(2*n1+1,2*n1+1)+v;