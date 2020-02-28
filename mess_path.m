function mypath = mess_path
%% Add all required directories to the MATLAB path
% Run this script to add all required functions and directories to the
% MATLAB path in order to rum M.E.S.S. functions and demos

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
mypath = genpath_exclude(pwd,{'.git','html'}); 
addpath(mypath);
if not(nargout)
    clear mypath;
end