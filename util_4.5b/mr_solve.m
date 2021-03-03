function [G, H]=mr_solve(Ap,A0,A_,B)

nd=length(Ap);
D = [zeros(nd,nd) Ap; eye(nd) zeros(nd,nd)];
E = [-A_ -A0; zeros(nd,nd) eye(nd)];

% [AA,BB,Q,Z] = qz(A,B) 
% Q*A*Z = AA, and Q*B*Z = BB
% [S, T, Q, Z] = qz(E,D,'complex');
[S, T, Q, Z] = qz(E,D,'real');
lam=diag(S)./diag(T);
% max(max(Q*E*Z-S))
% max(max(Q*D*Z-T))
% max(max(Q*E-S*Z'))
% max(max(Q*D-T*Z'))

SELECT=(abs(eig(S,T))<=1.000001);%(abs(lam)<=1.000001);
[SS,TS,QS,ZS] = ordqz(S,T,Q,Z,SELECT);
% [SS,TS,QS,ZS] = ordqz(S,T,Q,Z,'udi');
ns=length(find(SELECT));
nf=length(D)-ns;

% max(max(QS*E-SS*ZS'))
% max(max(QS*D-TS*ZS'))
ZZ=ZS';
Z11=ZZ(1:ns,1:ns);
Z12=ZZ(1:ns,ns+1:end);
Z21=ZZ(ns+1:end,1:ns);
Z22=ZZ(ns+1:end,ns+1:end);
G = -inv(Z22)*Z21;
H = -inv(Ap*G+A0)*B;   
