function eqn = getrail(k)

%% function eqn = getrail(k) reads in data from rail example
% s.a. Oberwolfach Model Reduction Benchmark Collection hosted at MORWiki
% (https://modelreduction.org)
%
% Input:
%   k       number of instance 
%           k=0 rail_371
%           k=1 rail_1357
%           k=2 rail_5177 
%           k=3 rail_20209
%           k=4 rail_79841
% Output:
% eqn       structure with members
%   E_      Sparse Matrix
%   A_      Sparse Matrix
%   B       Dense Matrix
%   C       Dense Matrix
%   haveE   indicates that E is not the identity

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
%

%% set path
switch k
    case 0
        example = 'rail_371';
    case 1
        example = 'rail_1357';
    case 2
        example = 'rail_5177';
    case 3  
        example = 'rail_20209';
    case 4
        example = 'rail_79841';
    otherwise
        error('MESS:error_arguments','k has to be 0, 1, 2, 3 or 4\n');
end

%% check path
path = fileparts(mfilename('fullpath'));

%% read matrices
eqn = load(strcat(path,'/',example),'E','A','B','C');
eqn.E_ = eqn.E;
eqn.A_ = eqn.A;
eqn = rmfield(eqn,'E');
eqn = rmfield(eqn,'A');
eqn.haveE = 1;



