function [eqn,opts,oper] = mul_A_pre_dae_2(eqn,opts,oper)
%  MUL_A_PRE_DAE_2 creates the hidden manifold projector used in
%  mul_A_dae_2.
%
% Input/Output: 
%  eqn, opts, oper  the usual equation, options and function handle
%                   structures; passed trough to enable storing
%                   helper data inside them. 
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
    
[eqn,opts,oper] = mul_Pi_pre(eqn,opts,oper);