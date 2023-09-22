function string = cleantex(input_string, backslashonly) %#ok<INUSD>
%% sanitizes input_string and escapes problematic characters with respect to TeX
%
%
%   input_string    String with potentially un-escaped characters
%
%   string          Escaped String for usage in TeX
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%
if nargin < 2
    input_string = strrep(input_string, '\n', '\newline{}');
    input_string = strrep(input_string, '_', '\_');
    input_string = strrep(input_string, '\', '\\');
    string = input_string;
else
    input_string = strrep(input_string, '\', '\\');
    string = input_string;
end
