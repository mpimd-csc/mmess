function [eqn, opts, oper] = mul_Pi_post(eqn,opts,oper)
% MUL_Pi multiplies with the hidden manifold projection matrix or it
% transpose. Note that the multiplication is actually implemented as the
% solution of a saddle point system.
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
if isfield(eqn,'Pcount')
    if eqn.Pcount
        eqn.Pcount = eqn.Pcount -1;
        if eqn.Pcount == 0
            eqn = rmfield(eqn, 'P_');
            eqn = rmfield(eqn, 'Pcount');
        end
    end
end