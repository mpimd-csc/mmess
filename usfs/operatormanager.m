function [oper, opts] = operatormanager(opts, name)
%% function [oper, opts] = operatormanager(opts, name)
%
%  Return structure with function handles that are implemented in folder
%  name. An error is thrown if the necessary function handles are not
%  given.
%
% Input
%   opts            the standard options structure required for logging
%   name            name of folder containing the function handle set
%
% Output
%   oper            struct, containing the function handles

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright (c) 2009-2023 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Check input parameter.
mess_assert(opts, ischar(name), ...
            'control_data', ...
            'input argument name has to be a string.');

%% Check path to function handles.
[fhpath, ~, ~] = fileparts(mfilename('fullpath'));
fhpath         = strcat(fhpath, filesep, name);

mess_assert(opts, any(exist(fhpath, 'dir')), ...
            'control_data', ...
            'there is no folder %s', fhpath);

%% Check for minimal function handle set.
funcs = { ...
         'mul_A', ...
         'mul_E', ...
         'size', ...
         'sol_A', ...
         'sol_E', ...
         'sol_ApE', ...
         'mul_ApE', ...
         'init', ...
         'init_res'};
for f = funcs
    mess_assert(opts, any(exist(strcat(fhpath, filesep, f{1}, '_', name, '.m'), ...
                                'file')), ...
                'check_data', ...
                'file %s does not exist', strcat(f{1}, '_', name, '.m'));
end

% Additional function check for state space transformed systems.
if any(strfind(name, 'state_space_transformed'))
    funcs = {'dss_to_ss', 'ss_to_dss'};
    for f = funcs
        mess_assert(opts, any(exist(strcat(fhpath, filesep, f{1}, '_', name, '.m'), ...
                                    'file')), ...
                    'check_data', ...
                    'file %s does not exist', strcat(f{1}, '_', name, '.m'));
    end
end

%% Create oper structure.
oper.name = name;

funcs = dir(fhpath);

for k = 3:length(funcs)
    % Sort out functions with wrong naming scheme and not M-files.
    [~, ~, file_ext] = fileparts(funcs(k).name);
    if not(any(strfind(funcs(k).name, name))) || not(strcmp(file_ext, '.m'))
        continue
    end

    fname = strrep(funcs(k).name, strcat('_', name, '.m'), '');

    % Put the existing function into the function handle set.
    eval(sprintf('oper.%s = @%s;', ...
                 fname, ...
                 strrep(funcs(k).name, '.m', '')));

    % Replace non-existing functions by do-nothing function.
    if not((any(strfind(fname, '_pre'))) || ...
           any(strfind(fname, '_post')) || ...
           exist(strcat(fhpath, filesep, fname, '_pre_', name, '.m'), ...
                 'file'))
        eval(sprintf('oper.%s = @mess_do_nothing;', ...
                     strcat(fname, '_pre')));
    end
    if not((any(strfind(fname, '_pre'))) || ...
           any(strfind(fname, '_post')) || ...
           exist(strcat(fhpath, filesep, fname, '_post_', name, '.m'), ...
                 'file'))
        eval(sprintf('oper.%s = @mess_do_nothing;', ...
                     strcat(fname, '_post')));
    end
end
