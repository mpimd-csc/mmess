function [ y, perm ] = mess_make_proper( x )
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
%          conjugation and with pairs as succesive entries
%  perm    the permutation applied to reorder the shifts, e.g. to
%          apply the correspunding reordering on the associated basis
%          vectors in IRKA as well.
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
%               2009-2020
%
if isempty(x), y=x; return; end

[m,n] = size(x);

if m < n, x = x'; end

idro = find(imag(x)==0.0);
idco = find(imag(x));

xr = x(idro);
xc = x(idco);

k = length(xc);
if rem(k,2)~=0, error('odd number of complex shifts detected'); end

% Sort real shifts
[xr,idr] = sort(xr);

% Sort complex shifts w.r.t. real part
[~,idcr] = sort(real(xc));
xc = xc(idcr); 
% Complex conjugated pairs ensured
xc(2:2:k,:) = conj(xc(1:2:k,:));    

% Sort complex shifts w.r.t. real part
[~, idcr2]=sort(real(xc));
xc = xc(idcr2);

y = [xr; xc];
perm = [idro(idr);idco(idcr(idcr2))];
end

