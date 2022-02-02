% mess_reset clears the work space variables, closes the open figures and
% clears the terminal.
%
% Note that only variable data is cleared, such that JIT compiler data can
% be preserved, which should speedup things over doing a manual 'clear all'.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2022 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

close all;

if exist('OCTAVE_VERSION', 'builtin')
    clear -variables;
else
    clearvars;
end

clc;