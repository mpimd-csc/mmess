function [Er,Ar,Br,Cr] = IRKA_FDM(n0,r,istest)
% IRKA_FDM computes a reduced order model via the IRKA method (see
% e.g. [1]) for a finite difference discretized convection
% diffusion model on the unit square described in [2]. 
%
% Usage: 
%    [Ar, Br, Cr] = IRKA_FDM(n0, r, test)
%
% Inputs
%
% n0          n0^2 gives the dimension of the original model, i.e. n0 is
%             the number of degrees of freedom, i.e. grid points, per
%             spatial direction 
%             (optional; defaults to 50)
% 
% r           desired dimension of the reduced order model 
%             (optional, defaults to 10)
%
% istest      flag to determine whether this demo runs as a CI test or 
%             interactive demo 
%             (optional, defaults to 0, i.e. interactive demo)
%
% Outputs
%
% Ar, Br, Cr  the reduced orde system matrices.
%
% References
% [1] A. C. Antoulas, Approximation of Large-Scale Dynamical Systems, Vol.
%     6 of Adv. Des. Control, SIAM Publications, Philadelphia, PA, 2005.
%     https://doi.org/10.1137/1.9780898718713.  
%
% [2] T. Penzl, Lyapack Users Guide, Tech. Rep. SFB393/00-33, 
%     Sonderforschungsbereich 393 Numerische Simulation auf massiv 
%     parallelen Rechnern, TU Chemnitz, 09107 Chemnitz, Germany, 
%     available from http://www.tu-chemnitz.de/sfb393/sfb00pr.html. (2000).

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

%%
if nargin<1
    n0=100;
end

E = speye(n0^2);
A = fdm_2d_matrix(n0,'10*x','100*y','0');
B = fdm_2d_vector(n0,'.1<x<=.3');
C = fdm_2d_vector(n0,'.7<x<=.9');C=C';
%%
if nargin<2
    opts.irka.r = 10; 
else
    opts.irka.r = r; 
end
if nargin < 3
    istest = 0;
end
opts.irka.maxiter =20;
opts.irka.shift_tol = 1e-2;
opts.irka.h2_tol = 1e-6;
if istest
    opts.irka.info = 1;
else
    opts.irka.info = 2;
end
opts.irka.init = 'logspace';

[Er,Ar,Br,Cr,~,~,~,~,~] = mess_tangential_irka(E,A,B,C,opts);

if istest
    [~, ID] = lastwarn;
    lastwarn('');
    if strcmp(ID,'MESS:IRKA:convergence')
        error('IRKA converged unexpectedly slow');
    end
end