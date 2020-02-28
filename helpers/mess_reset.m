% mess_reset clears the work space variables, closes the open figures and
% clears the terminal.
%
% Note that only variable data is cleared, such that JIT compiler data can
% be preserved, which should speedup things over doing a manual 'clear all'.

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
close all;

if exist('OCTAVE_VERSION', 'builtin')
    clear -variables;
else
    clearvars;
end

clc;