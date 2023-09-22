function mypath = mess_path(prototypes, force)
%% Add all required directories to the MATLAB path
% Run this script to add all required functions and directories to the
% MATLAB path in order to run M.E.S.S. functions and demos or
% generate a list of directories for permanent addition to your
% MATLAB path.
%
% Calls:
%   mess_path
%   pathlist = mess_path
%
% on the development version you may want to run
%
%   mess_path(true)
%   pathlist = mess_path(true)
%
% to also add the _prototypes folder.
%
% If you want to force-add the path (e.g. because you want to add
% the _prototypes folder after having added the base folder without
% it you may run
%
%   mess_path(true, true)
%   pathlist = mess_path(true, true)
%
% Note that we only add folders to the path, i.e. something like
%
%   mess_path(false, true)
%
% will NOT remove the _prototypes folders if they have been added
% before.
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Let's check if another version is already on the path
if exist('mess_lyap', 'file')
    if (nargin > 1) && force
        err_warn = @warning;
    else
        err_warn = @error;
    end

    err_warn('MESS:path_exists', ...
             ['It seems like M-M.E.S.S. is already ' ...
              'on your MATLAB search path. You should remove ' ...
              'the existing instance from the path to ' ...
              'avoid version conflicts.']);
end

%% Now generate and add this versions path
if (nargin > 0) && prototypes
    mypath = genpath_exclude(pwd, ...
                             {'.git', 'html', '_packages', 'resources'});
else
    mypath = genpath_exclude(pwd, ...
                             {'.git', 'html', '_prototypes', '_packages', 'resources'});
end

addpath(mypath);

if not(nargout)
    clear mypath;
end
