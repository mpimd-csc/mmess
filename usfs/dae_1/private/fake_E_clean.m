function eqn = fake_E_clean(eqn)
if eqn.Ecount > 1
    eqn.Ecount = eqn.Ecount -1;
else
    eqn = rmfield(eqn,{'E_','Ecount'});
end