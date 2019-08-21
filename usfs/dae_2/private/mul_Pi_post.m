function [eqn, opts, oper] = mul_Pi_post(eqn,opts,oper)
% MUL_Pi multiplies with the hidden manifold projection matrix or it
% transpose. Note that the multiplication is actually implemented as the
% solution of a saddle point system.
% 

if isfield(eqn,'Pcount')
    if eqn.Pcount
        eqn.Pcount = eqn.Pcount -1;
        if eqn.Pcount == 0
            eqn = rmfield(eqn, 'P_');
            eqn = rmfield(eqn, 'Pcount');
        end
    end
end