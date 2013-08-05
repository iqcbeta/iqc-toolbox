function [w,X,x]=iqc_delay1(v,T0,a)
%  function [w,X,x]=iqc_delay1(v,T0,a)
%
% Defines IQCs for the relation w=e^{-sT}v,
% where 0<T<T0 is an uncertain time delay.
%
% The IQCs have the form 
%
% (H*v)'*X*(H*v)>(w-v)'*X*(w-v)
%
% where H(s)=2sT0(sT0+srt(12.5))/((sT0)^2+asT0+b),
%       b=sqrt(50), a=sqrt(2b+6.5)
%
% and X=x0+x1/(s+a(1))+....+xN/(s+a(N))>0
%     a(i)>0, i=1,..,N
%     xi arbitrary scalar variable 
%
% Default a=[]
%         T0=1
%
% Work by UlfJ June 12, 1997


if nargin<1
   error('input must be specified')
end
if size(v,1)~=1, error('only scalar signals are allowed'), end
if nargin<2, T0=1; end
if nargin<3, a=[]; end

s=tf([1 0],1);
w=signal;
x{1}=rectangular;
X=x{1};
N=length(a);
for k=1:N
  x{k+1}=rectangular;
  X=X+x{k+1}*(1/(s+a(k)));
end
X>0;
b=sqrt(50);
a=sqrt(2*b+6.5);
H=2*T0*s*(T0*s+sqrt(12.5))/(T0^2*s*s+a*T0*s+b);
(H*v)'*X*(H*v)>(w-v)'*X*(w-v);

