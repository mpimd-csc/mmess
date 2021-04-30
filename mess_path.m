function mypath = mess_path(prototypes)
%% Add all required directories to the MATLAB path
% Run this script to add all required functions and directories to the
% MATLAB path in order to run M.E.S.S. functions and demos or
% generate a list of directories for permanent addition to your
% MATLAB path.
%
% Calls:
%   messpath
%   pathlist = mess_path
%
% on the development version you may want to run
%
%   messpath(true)
%   pathlist = mess_path(true)
%
% to also add the prototypes folder.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if (nargin > 0) && prototypes
    mypath = genpath_exclude(pwd,{'.git','html'});
else
    mypath = genpath_exclude(pwd,{'.git','html','_prototypes','_packages'});
end

addpath(mypath);

if not(nargout)
    clear mypath;
end
