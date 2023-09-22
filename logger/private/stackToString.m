function StackString = stackToString(stack)
% creates a compact backtrace from a dbstack

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
StackString = '';
for i = 2:length(stack)
    StackString = strcat(StackString, '>>', stack(i).file, '(', int2str(stack(i).line), ')');
end
