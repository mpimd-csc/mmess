function out = mess_string(in)
% function out = mess_string(in)
%
% converts input into a string, or array of strings if it was an array or
% cell itself
%
% mess_string falls back to string if that exists and passed input strings
% directly to the output.
%
% This is mostly intended to add the string function to octave, where it
% does not exist (yet).
%

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

out = '';
if isa(in, 'char')
    out = in;
elseif exist('string', 'file')
    out = string(in);
elseif isa(in, 'cell')
    for sub = in
        out = [out, mess_string(sub)]; %#ok<AGROW>
    end
elseif isa(in, 'double') || isa(in, 'single')
    for i = 1:length(in)
        out = sprintf('%s, %g', out, in);
    end
elseif isa(in, 'int8') || isa(in, 'int16') || isa(in, 'int32') || ...
        isa(in, 'int64') || isa(in, 'uint8') || isa(in, 'uint16') || ...
        isa(in, 'uint32') || isa(in, 'uint64')
    for i = 1:length(in)
        out = sprintf('%s, %d', in);
    end
elseif isa(in, 'logical')
    Logicals = {'False', 'True'};
    out = Logicals{in(1) + 1};
    for i = 2:length(in)
        out = [out Logicals(in(i) + 1)]; %#ok<AGROW>
    end
end
