function [rs,ru]=separate(g)
% function r=separate(g)
% given a scalar system g without poles on jR,
% (and without repeated poles)
% find a stable system rs
% and an antistable sustem ru such that g=rs+ru
% (the direct terms of rs and ru must be equal)

[A,B,C,D]=ssdata(g);
d=D/2;
[v,l]=eig(A);
l=diag(l);
n=length(l);
ns=length(find(real(l)<0));
v=[v(:,find(real(l)<0)) v(:,find(real(l)>0))];
vi=inv(v);
c=C*v;
a=vi*A*v;
b=vi*B;
[nrs,drs]=ss2tf(a(1:ns,1:ns),b(1:ns),c(1:ns),d);
rs=tf(nrs,drs);
[nru,dru]=ss2tf(a(ns+1:n,ns+1:n),b(ns+1:n),c(ns+1:n),d);
ru=tf(nru,dru);
