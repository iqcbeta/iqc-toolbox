function [w,X,x]=iqc_window(v,T0,a)
% function [w,X,x]=iqc_window(v,T0,a)
%
% Defines IQCs for the relation
%
%   w=((e^{-sT}-1)/s)*v,
%
% where 0<T<T0 is an uncertain time delay. This corresponds to
% convolution with a rectangular window in the time domain.
%
% The IQCs have the form
%
% (H*v)'*X*(H*v)>w'*X*w
%
% where H=2T0(s*T0+sqrt(12.5))/((sT0)^2+asT0+b),
%       b=sqrt(50), a=sqrt(2b+6.5)
%
% and X(s)=x0+x1/(s+a(1))+....+xN/(s+a(N))>0
%     a(i)>0, i=1,..,N
%     xi arbitrary scalar variable
%
% Default a=[]
%         T0=1
%
% Outputs:   [x0]
%          x=[ :]
%            [xN]
%          X(s)
%
% Work by UlfJ Oct 21, 1997
% Last modified by cmj on 2013/5/5

if nargin<1
    disp_str(1)
end
if size(v,1)~=1
    disp_str(8,'scalar')
end
if nargin<2
    T0=1;
end
if nargin<3
    a=[];
end

global ABST

switch ABST.systemtype
    case 'continuous'
        s=tf([1 0],1);
        w=signal;
        x{1}=rectangular;
        X=x{1};
        N=length(a);
        for k=1:length(a)
            x{k+1}=rectangular;
            X=X+x{k+1}*(1/(s+a(k)));
        end
        X>0;
        b=sqrt(50);
        a=sqrt(2*b+6.5);
        H=2*T0*(T0*s+sqrt(12.5))/(T0^2*s*s+a*T0*s+b);
        (H*v)'*X*(H*v)>w'*X*w;
    case 'discrete'
        disp_str(70,'iqc_window','discrete')
end
