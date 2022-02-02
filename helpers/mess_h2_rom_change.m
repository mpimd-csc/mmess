function [romchg] = mess_h2_rom_change(E1,A1,B1,C1,E2,A2,B2,C2,rel)
% [romchg] = mess_h2_rom_change(E1,A1,B1,C1,E2,A2,B2,C2,rel)
%
% computes the (relative) difference of two stable systems in the H2 norm.
%
% Inputs:
% E1,A1,B1,C1,E2,A2,B2,C2  The system matrices (E1,E2 invertible)
% rel                      indicator whether the relative or absolute norm
%                          is desired.
%
% Output:
% romchg                   the computed H2-norm difference
%

%
% This file is part of the M-M.E.S.S. project 
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if nargin<8, error('to few inputs'); end
if nargin<9, rel=0; end

E = blkdiag(E1,E2);
A = blkdiag(A1,A2);
B = [B1; B2];
C = [C1, -C2];

if exist('lyap','file')
    X = lyap(A,B*B',[],E);

else
    B = E\B;
    X = lyap2solve(E\A,B*B');
end
nrm = sqrt(trace(C*(X*C')));
if rel
    if exist('lyap','file')
        X1 = lyap(A1,B1*B1',[],E1);
    else
        B1 = E1\B1;
        X1 = lyap2solve(E1\A1,B1*B1');
    end
    nrm1 = sqrt(trace(C1*(X1*C1')));
else
    nrm1 = 1.0;
end
romchg = nrm/nrm1;
