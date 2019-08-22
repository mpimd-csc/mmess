function [Er,Ar,Br,Cr] = IRKA_rail(k,r,istest)
% Computes a loclly H2-optimal reduced order model of order r for
% the selective cooling of Steel profiles application described in
% [1,2,3] via the tangential IRKA method.
%
% Usage: IRKA_rail(k,r,istest)
%
% Inputs: 
% 
% k           refinement level of the model to use 
%             (1-4, i.e. 1357-79841Dofs)
%             (optinal, defaults to 1)
%
% r           desired dimension of the reduced order model 
%             (optional, defaults to 10)
%
% istest      flag to determine whether this demo runs as a CI test or 
%             interactive demo
%             (optional, defaults to 0, i.e. interactive demo)
%
% References:
% [1] J. Saak, Effiziente numerische Lösung eines
%     Optimalsteuerungsproblems für die Abkühlung von Stahlprofilen,
%     Diplomarbeit, Fachbereich 3/Mathematik und Informatik, Universität
%     Bremen, D-28334 Bremen (Sep. 2003).   
%
% [2] P. Benner, J. Saak, A semi-discretized heat transfer model for
%     optimal cooling of steel profiles, in: P. Benner, V. Mehrmann, D.
%     Sorensen (Eds.), Dimension Reduction of Large-Scale Systems, Vol. 45
%     of Lect. Notes Comput. Sci. Eng., Springer-Verlag, Berlin/Heidelberg,
%     Germany, 2005, pp. 353–356. https://doi.org/10.1007/3-540-27909-1_19. 
%
% [3] J. Saak, Efficient numerical solution of large scale algebraic matrix
%     equations in PDE control and model order reduction, Dissertation,
%     Technische Universität Chemnitz, Chemnitz, Germany (Jul. 2009).  
%     URL http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642
%
% [4] S. Gugercin, A. C. Antoulas, C. Beattie, H2 model reduction
%     for large-scale linear dynamical systems, SIAM J. Matrix
%     Anal. Appl. 30 (2) (2008) pp. 609–638. 
%     https://doi.org/10.1137/060666123. 

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
    k=1;
end
eqn = getrail(k);
%%
if nargin<2
    opts.irka.r = 20; 
else
    opts.irka.r = r; 
end

if nargin < 3
    istest = 0;
end
opts.irka.maxiter =20;
opts.irka.shift_tol = 1e-3;
opts.irka.h2_tol = 1e-6;
if istest
    opts.irka.info = 1;
else
    opts.irka.info = 2;
end
opts.irka.init = 'logspace';

[Er,Ar,Br,Cr] = mess_tangential_irka(eqn.E_,eqn.A_,eqn.B,eqn.C,opts);

if istest
    [~, ID] = lastwarn;
    lastwarn('');
    if strcmp(ID,'MESS:IRKA:convergence')
        error('IRKA converged unexpectedly slow');
    end
end