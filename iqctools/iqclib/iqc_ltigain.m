function [w,p,q,x,y,z]=iqc_gain(v,a)
% function [w,p,q,x,y,z]=iqc_ltigain(v,a)
%
% defines iqc's for the relation 
% w(t)=delta*v(t), 
% where delta is a constant scalar
% between -1 and 1
% (v,w could be vectors)
%
% the iqc's have the form
% v'*p*v>w'*p*w
% where p=x0+x1/(s+a(1))+x2/(s+a(2))+...>0,
% and also 
% v'*q*w>0
% where 
% q/2=z0+(y1*s+a(1)*z1)/(s^2-a(1)^2)+(y2*s+a(2)*z2)/(s^2-a(2)^2)+...
% xi arbitrary square matrices,
% yi arbitrary symmetric,
% zi arbitrary skew symmetric
%
% default a=1
%
% Written by ameg@mit.edu,  last modified October 13, 1997
if nargin<2, a=1; end

m=size(v,1);  
n=length(a);

s=tf([1 0],1);

w=signal(m);

if m>1,
   x0=rectangular(m);
   y0=symmetric(m);
   z0=skew(m);
   p=x0;
   q=z0;
   v'*(x0*v)>w'*x0*w;
   v'*z0*w==0;
   for i=1:n,
      x{i}=rectangular(m);
      y{i}=symmetric(m);
      z{i}=skew(m);
      Ga=abst(1/(s+a(i)));
      va=Ga*v;
      wa=Ga*w;
      p=p+x{i}*Ga;
      q=q+y{i}*Ga-Ga'*y{i}+z{i}*Ga+Ga'*z{i};
      v'*(x{i}*va)>w'*x{i}*wa;
      v'*y{i}*wa>va'*(y{i}*w);
      v'*z{i}*wa+va'*z{i}*w>0;
   end
else
   x0=rectangular;
   y0=symmetric;
   p=x0;
   q=0;
   v'*(x0*v)>w'*x0*w;
   for i=1:n,
      x{i}=rectangular;
      y{i}=symmetric;
      Ga=abst(1/(s+a(i)));
      va=Ga*v;
      wa=Ga*w;
      p=p+x{i}*Ga;
      q=q+y{i}*Ga-Ga'*y{i};
      v'*(x{i}*va)>w'*x{i}*wa;
      v'*y{i}*wa==va'*(y{i}*w);
   end
end
p>0;