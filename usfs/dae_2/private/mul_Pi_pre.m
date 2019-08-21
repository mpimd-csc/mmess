function [eqn, opts, oper] = mul_Pi_pre(eqn,opts,oper)
% MUL_Pi multiplies with the hidden manifold projection matrix or it
% transpose. Note that the multiplication is actually implemented as the
% solution of a saddle point system.
% 

if not(isfield(eqn,'P_'))
    eqn.P_ = eqn.A_;
    eqn.P_(1:eqn.st, 1:eqn.st) = eqn.E_(1:eqn.st,1:eqn.st);
    eqn.Pcount = 1;
else
    eqn.Pcount = eqn.Pcount +1;
end

