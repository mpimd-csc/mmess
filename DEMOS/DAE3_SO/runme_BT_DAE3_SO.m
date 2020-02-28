%  Driver script for the results presented in:
%
%  J. Saak, M. Voigt, Model reduction of constrained mechanical systems
%  in M-M.E.S.S., IFAC-PapersOnLine 9th Vienna International Conference
%  on Mathematical Modelling MATHMOD 2018, Vienna, Austria, 21–23
%  February 2018 51 (2) (2018) 661–666.
%  https://doi.org/10.1016/j.ifacol.2018.03.112.      

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
model={{'Stykel_large', 1e-4, 200, 500}, ...
       {'Truhar_Veselic', 1e-4, 250, 300}, ...
       {'TV2', 1e-3, 300, 300},...
       {'TV2', 1e-3, 300, 500}...
      };
l = length(model);
for i=1:l
    BT_DAE3_SO(model{i}{1},model{i}{2},model{i}{3},model{i}{4})
    fprintf('\nPress any key to continue\n');
    if (i<l), pause; end
end
