function mess_websave(filename, downloadurl)
%%% mess_websave is a websave wrapper that falls back to urlwrite
%   on older MATLAB and OCTAVE where websave is not available.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
if exist('websave', 'file')
    websave(filename, downloadurl);
else
    urlwrite(downloadurl, filename); %#ok<URLWR>
end
