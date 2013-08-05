function z=lti_a(num,den,a,x,v)
% function z=lti_a(num,den,a,x,v)
% this function aids efficient iqc-programming 
% for some uncertain LTI systems:
% without adding new states, it forms signal
% z=(1/(s+a))*(num(s)/den(s))*u,
% given x: [y,x]=lti_full(num,den,u)
%       v:     v=(1/(s+a))*u
% this works for VECTOR u 
%
% Written by ameg@mit.edu, last modified Dec. 07, 1997

n=length(den)-1;
g=polyval(num,-a)/polyval(den,-a);
num=[zeros(1,n+1-length(num)) num];
r=deconv(num-g*den,[1 a]);
z=g*v;
nr=length(r);
for k=1:nr,
   z=z+r(nr-k+1)*x{k};
end
