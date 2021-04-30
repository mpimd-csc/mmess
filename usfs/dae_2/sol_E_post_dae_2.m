function [ eqn, opts, oper ] = sol_E_post_dae_2( eqn, opts, oper )
%% function post finalizes data and/or functions
%
% Input:
%    eqn    struct contains data for equations
%
%    opts   struct contains parameters for the algorithm
%
%    oper   struct contains function handles for operation with A
%
% Output:
% eqn
% opts
% oper

%
% This file is part of the M-M.E.S.S. project
% (http://www.mpi-magdeburg.mpg.de/projects/mess).
% Copyright © 2009-2021 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

  if(not(isfield(eqn, 'Scount'))) || not(isnumeric(eqn.Scount))
      error('MESS:error_arguments', ['field eqn.Scount is not defined. Did ' ...
                        'you forget to run mul_E_pre?']);
  end
  if eqn.Scount>1
    eqn.Scount=eqn.Scount-1;
  else
    eqn=rmfield(eqn,'S_');
    eqn=rmfield(eqn,'Scount');
  end
end
