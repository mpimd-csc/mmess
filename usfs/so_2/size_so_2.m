function n = size_so_2(eqn, opts)%#ok<INUSD>

% function n = size_so_2(eqn, opts)
%
% Call help mess_usfs_so_2 to see the description of the second order
% system and its transformed first order system
%
%
% This function returns the number of rows of the matrices A and E.
%
%    Inputs:
%
%    eqn       structure containing field 'K_'
%    opts      structure containing parameters for the algorithm
%
%    Output:
%
%    n         double number of rows of matrix V_ in structure eqn
%
% This function does not use other so1 functions.

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% check data
if not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_))
    error('MESS:error_arguments',...
    'Field eqn.K_ is not defined or corrupted');
end

%% compute size
n = 2*size(eqn.K_,1);

end
