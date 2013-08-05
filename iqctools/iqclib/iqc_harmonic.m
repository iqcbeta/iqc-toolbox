function [w,p,q,x]=iqc_gain(v,w0,a)
% function [w,p,q,x]=iqc_ltigain(v,w0,a)
%
% defines iqc's for the relation 
% w(t)=delta(t)*v(t), 
% where delta(t)=cos(w0*t+p0)
% (v,w could be vectors)
%
% the iqc's have the form
% v'*q*v>w'*p*w
% where 
%   p=x0+x1/(s+a(1))+x2/(s+a(2))+...>0,
%   q=x0+x1(s+a{1})/(s^2+2*a{1}*s+a{1}^2+w0^2)+...
% xi arbitrary square matrices,
% default w0=1, a=1
%
% Written by ameg@mit.edu,  last modified November 29, 1997
if nargin<2, w0=1; end
if nargin<3, a=1; end

m=size(v,1);  
n=length(a);
s=tf([1 0],1);
w=signal(m);
x0=rectangular(m);
p=x0;
q=x0;
v'*(x0*v)>w'*x0*w;
for i=1:n,
   x{i}=rectangular(m);
   Ga=abst(1/(s+a(i)));
   Gaw=abst((s+a(i))/(s*s+2*a(i)*s+a(i)^2+w0^2));
   p=p+x{i}*Ga;
   q=q+x{i}*Gaw;
   v'*x{i}*(Gaw*v)>w'*x{i}*(Ga*w);
end
p>0;