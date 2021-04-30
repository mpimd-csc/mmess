function [ eqn, opts, oper ] = sol_E_pre_dae_2( eqn, opts, oper )
%% function pre initializes data and/or functions
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

  alpha=-1/50;
  if isfield(eqn,'st')&&isnumeric(eqn.st)
    st=eqn.st;
  else
    error('MESS:wrong_arguments','missing or corrupted field st detected');
  end
  if not(isfield(eqn,'S_'))
    if(not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_))...
            || not(isfield(eqn,'A_'))) || not(isnumeric(eqn.A_))
        error('MESS:error_arguments','field eqn.E_ or eqn.A_ is not defined or corrupted');
    end
    eqn.S_=alpha*eqn.A_;
    eqn.S_(1:st,1:st)=eqn.E_(1:st,1:st);
    eqn.Scount=1;
  else
    if(not(isfield(eqn, 'Scount'))) || not(isnumeric(eqn.Scount))
        error('MESS:error_arguments', ['field eqn.Scount is not defined. Did ' ...
                        'you forget to run sol_E_pre?']);
    end
    eqn.Scount=eqn.Scount+1;
  end
end

