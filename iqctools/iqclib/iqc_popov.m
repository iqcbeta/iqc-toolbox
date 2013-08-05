function [w,x,y]=iqc_popov(v,alpha,beta)
% function [w,x,y]=iqc_popov(v,alpha,beta)
%
% Popov IQC for the sector nonlinearity
% w(t)=f(v(t)), where f is in the [alpba,beta] sector:
%
% (w-a*v)'*x*(b*v-w)+w'*y*(dv/dt)>0, where x>0
%
% default alpha=0, beta=1
%
% Written by ameg@mit.edu,  last modified October 13, 1997
if nargin<2, 
   a=0;
   b=1;
elseif nargin<3
  if length(alpha)<2
     a=alpha;
     b=1;
  else
     a=alpha(1);
     b=alpha(2);
  end
else
  a=alpha;
  b=beta; 
end

if ~(a<b), error('empty sector'); end

u=derivative(v);
w=signal;
x=symmetric;
y=symmetric;

x>0;
(w-a*v)'*x*(b*v-w)+w'*y*u>0;