function eqn = fake_E(eqn)
if not(isfield(eqn,'Ecount'))
    eqn.Ecount = 1;
    % E = [ I 0 ]
    %     [ 0 0 ]
    eqn.E_=sparse(1:st,1:st,ones(st, 1),n,n,st);
else
    eqn.Ecount = eqn.Ecount + 1;
end
