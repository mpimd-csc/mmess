function [eqn,opts,oper] = mul_A_post_dae_2(eqn,opts,oper)
    
[eqn,opts,oper] = mul_Pi_post(eqn,opts,oper);