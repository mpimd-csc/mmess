function [y,yr] = impeuler(E,A,B,C,Er,Ar,Br,Cr,tau,tmin,tmax,x0,xr0,alpha)
% Simple implicit Euler implementation for validation of the DAE2
% MESS open loop example via a basic step response computation
%
%  [y,yr] = impeuler(E,A,B,C,Er,Ar,Br,Cr,tau,tmin,tmax,x0,xr0,alpha)
%
% INPUTS:
% E,A,B,C      The original system matrices
% Er,Ar,Br,Cr  The reduced order system matrices
% tau          time step size
% tmin         start time
% tmax         end time
% x0, x0r      initial states of the full and reduced order models
% alpha        step height for the input step function
%
% OUTPUT:
% y,yr         outputs of the full and reduced systems in [tmin,tmax]
%


%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%


  [L,U,P,Q] = lu(E-tau*A);
  [Lr,Ur,Pr] = lu(Er-tau*Ar);

 ntau=ceil((tmax-tmin)/tau);
 y=zeros(size(C,1),ntau);
 yr=y;
 for i=1:ntau
    if not(mod(i,ceil(ntau/10)))
      fprintf('\r Implicit Euler step %d / %d',i,ntau);
    end
    if i<(0.1*ntau)
        x=Q*(U\(L\(P*(E*x0))));
        xr=Ur\(Lr\(Pr*(Er*xr0)));
    else
        salpha=smoother(alpha,i,ntau);
        x=Q*(U\(L\(P*(E*x0+(salpha*tau)*sum(B,2)))));
        xr=Ur\(Lr\(Pr*(Er*xr0+(salpha*tau)*sum(Br,2))));
    end
    y(:,i)=C*x;
    yr(:,i)=Cr*xr;
    x0=x;
    xr0=xr;
 end
  fprintf('\n\n');
end

