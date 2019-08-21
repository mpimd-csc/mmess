function [TL,TR] = square_root_method_SO(M, max_ord, tol, info, U,V)
[U0,S0,V0] = svd(V'*M*U,0);


s0=diag(S0);
ks=length(s0);
k0=1;
while (k0<=ks) && (s0(k0)>s0(1)*tol)
    k0=k0+1;
end


r= min([max_ord k0]);
if info>0
    fprintf(1,'reduced system order: %d  (max possible/allowed: %d/%d)\n\n',r,ks,max_ord);
end

sigma_r=diag(S0(1:r,1:r));

VB = U*V0(:,1:r);
VC = V*U0(:,1:r);

TR = VB*diag(ones(r,1)./sqrt(sigma_r));
TL = VC*diag(ones(r,1)./sqrt(sigma_r));
end