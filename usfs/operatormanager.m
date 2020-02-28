function oper = operatormanager(name)
%% function oper = operatormanager(name)
%
%  Return structure with functionhandles that are implemented in folder
%  name. An error is thrown if the necessary function handles are not
%  given.
%
% Input
%   name            name of folder containing the function handle set
%
% Output
%	oper            struct, containing the function handles

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

%% Check input parameter.
assert(ischar(name), ...
    'MESS:control_data', ...
    'input argument name has to be a string.');

%% Check path to function handles.
[fhpath, ~, ~] = fileparts(mfilename('fullpath'));
fhpath         = strcat(fhpath, '/', name);

assert(any(exist(fhpath, 'dir')), ...
    'MESS:control_data', ...
    'theres no folder %s', fhpath);

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
    assert(any(exist(strcat(fhpath, '/', f{1}, '_', name, '.m'), ...
        'file')), ...
        'MESS:check_data', ...
        'file %s does not exist', strcat(f{1}, '_', name, '.m'));
end

% Additional function check for state space transformed systems.
if any(strfind(name, 'state_space_transformed'))
    funcs = {'dss_to_ss', 'ss_to_dss'};
    for f = funcs
        assert(any(exist(strcat(fhpath, '/', f{1}, '_', name, '.m'), ...
            'file')), ...
            'MESS:check_data', ...
            'file %s does not exist', strcat(f{1}, '_', name, '.m'));
    end
end

%% Create oper structure.
oper.name = name;

funcs = dir(fhpath);

for k = 3:length(funcs)
    % Sort out functions with wrong naming scheme and not M-files.
    [~,~,file_ext] = fileparts(funcs(k).name);
    if not(any(strfind(funcs(k).name, name))) || not(strcmp(file_ext,'.m'))
        continue;
    end
    
    fname = strrep(funcs(k).name, strcat('_', name, '.m'), '');
    
    % Put the existing function into the function handle set.
    eval(sprintf('oper.%s = @%s;', ...
        fname, ...
        strrep(funcs(k).name, '.m', '')));
    
    % Replace non-existing functions by do-nothing function.
    if not((any(strfind(fname, '_pre'))) ...
            || any(strfind(fname, '_post')) ...
            || exist(strcat(fhpath, '/',fname, '_pre_', name, '.m'), ...
            'file'))
        eval(sprintf('oper.%s = @mess_do_nothing;', ...
            strcat(fname, '_pre')));
    end
    if not((any(strfind(fname, '_pre'))) ...
            || any(strfind(fname, '_post')) ...
            || exist(strcat(fhpath, '/',fname, '_post_', name, '.m'), ...
            'file'))
        eval(sprintf('oper.%s = @mess_do_nothing;', ...
            strcat(fname, '_post')));
    end
end
