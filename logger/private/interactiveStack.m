function StackString = interactiveStack(stack)
% creates a backtrace from a dbstack with opentoline ref links of the
% backtrace

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

if length(stack) <= 2
    arrow = '>';
else
    arrow = '\x21b3';
end
if exist('OCTAVE_VERSION', 'builtin')
    % octave does not support opentoline
    for i = 2:length(stack)
        StackString{i - 1} = sprintf('%s In %s (line %s)\n', ...
                                     '>', stack(i).name, ...
                                     int2str(stack(i).line)); %#ok<*AGROW>
    end
else
    for i = 2:length(stack)
        StackString{i - 1} = sprintf('%s[\b In <a href="matlab:opentoline(''%s'',%s)">%s (line %s)</a>]\b\n', ...
                                     arrow, cleantex(which(stack(i).file), true),  int2str(stack(i).line), ...
                                     stack(i).name, int2str(stack(i).line));
    end
end
