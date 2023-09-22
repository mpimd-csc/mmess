function [y, perm] = mess_make_proper(x)
% MESS_MAKE_PROPER ensures the input vector to be a proper set of shifts
%
%   [ y, perm ] = mess_make_proper( x )
%
%  The function checks if the vector contains a set of shifts closed under
%  complex conjugation and all complex pairs are subsequent entries.
%
%  Input:
%
%  x       a vector of potentially complex shift parameters (e.g. for ADI, or IRKA)
%
%  Output:
%
%  y       a vector of proper shifts, i.e. closed under complex
%          conjugation and with pairs as successive entries
%  perm    the permutation applied to reorder the shifts, e.g. to
%          apply the corresponding reordering on the associated basis
%          vectors in IRKA as well.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
opts = struct;
if isempty(x)
    y = x;
    return
end

[m, n] = size(x);

if m < n
    x = x';
end

idro = find(imag(x) == 0.0);
idco = find(imag(x));

xr = x(idro);
xc = x(idco);

k = length(xc);
if not(rem(k, 2) == 0)
    mess_err(opts, 'error_arguments', 'odd number of complex shifts detected');
end

% Sort real shifts
[xr, idr] = sort(xr);

% Sort complex shifts w.r.t. real part
[~, idcr] = sort(real(xc));
xc = xc(idcr);
% Complex conjugated pairs ensured
xc(2:2:k, :) = conj(xc(1:2:k, :));

% Sort complex shifts w.r.t. real part
[~, idcr2] = sort(real(xc));
xc = xc(idcr2);

y = [xr; xc];
perm = [idro(idr); idco(idcr(idcr2))];
end
