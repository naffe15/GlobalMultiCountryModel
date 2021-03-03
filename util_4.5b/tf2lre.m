function [NUMNUM, DENDEN] = tf2lre(a, b, C, b1),

if nargin<4,
  b1=b;
end
na = length(a);
ns = length(a)-1;
nu = length(b);
nf = length(C)-na;
nd = ns+nu+nf+1;

G0=eye(nd);
for j=1:nf+1,
  G0(j,j:na+nu+j-1)=[a -b];
  G0(j,na+j:na+nf )=0;
end

% G0 = [1     a(2)  a(3)  a(4)  0     0     0     0     0     0    ;
%       0     1     a(2)  a(3)  a(4)  0     0     -b(3) 0     0    ;
%       0     0     1     a(2)  a(3)  a(4)  0     -b(2) -b(3) 0    ;
%       0     0     0     1     a(2)  a(3)  a(4)  -b(1) -b(2) -b(3);
%       zeros(3,4)              eye(3)            zeros(3,3);
%       0     0     0     0     0     0     0     1    0     0    ;
%       zeros(2,8)                                eye(2)];
    
G_ = [zeros(nf+1,nd);
  zeros(ns,nf) -eye(ns) zeros(ns,nu+1);
  zeros(1,nd);
  zeros(nu-1,nf+na) -eye(nu-1) zeros(nu-1,1)];

G = -inv(G0)*G_;
    
H0 = zeros(nd,1);
H0(nf+na+1) = -1;

H= -inv(G0)*H0;

A_ = G_;
Ap = [zeros(nf,1) -eye(nf) zeros(nf,ns+nu);
  zeros(na+nu,nd)];

B = H0;

% C=poly([10 4 1.5 0.95 0.6 0.5]);
% % C=poly([9 4.5 1.35 0.95 0.6 0.5]);
% C=C./C(4);

A0 = eye(nd);
A0(nf+1,:)=[C -b1];
% A0 = [1     0     0     0     0     0     0     0     0     0    ;
%       0     1     0     0     0     0     0     0     0     0    ;
%       0     0     1     0     0     0     0     0     0     0    ;
%       C                                         -b1 -b2 -b3;
%       zeros(3,4)              eye(3)            zeros(3,3);
%       0     0     0     0     0     0     0     1    0     0    ;
%       zeros(2,8)                                eye(2)];

% [H -inv(Ap*G+A0)*B ]   

D = [zeros(nd,nd) Ap; eye(nd) zeros(nd,nd)];
E = [-A_ -A0; zeros(nd,nd) eye(nd)];

% [V,LAM] = eig(E,D);
% max(max(D*V- LAM*E*V))

[Gok,Hok]=mr_solve(Ap,A0,A_,B);
max(max(Gok-G));

if matlab_ver_less_than('7.0')
  [x, fmin, EXITFLAG] = fminsearch('mr_sid_fn',b1,[], G, H, Ap, A0, A_, B);
else
  [x, fmin, EXITFLAG] = fminsearch(@(x) mr_sid_fn(x, G, H, Ap, A0, A_, B), b1);
end
if EXITFLAG==1,
  disp(['Minimum norm of the deviation of G is: ',num2str(exp(fmin))])
else
  disp(['Optimization did not converge'])
  disp(['Current norm of the deviation of G is: ',num2str(exp(fmin))])
end

DENDEN = C;
NUMNUM = x';

