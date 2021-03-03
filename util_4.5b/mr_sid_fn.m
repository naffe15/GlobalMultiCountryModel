function [ln, Gx, Hx] = mr_sid_fn(x, G, H, Ap, A0, A_, B),

nf=rank(Ap);
Ax=A0;
Ax(nf+1,end-length(x)+1:end)=x';

[Gx,Hx]=mr_solve(Ap,Ax,A_,B);

try
nn=norm([G H]-[Gx Hx]);
catch
  nn=Inf;
end
ln=log(nn);
