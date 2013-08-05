function [w,K1,M1,K2,M2]=iqc_slowtv(v,d,a,D)
% 
% function [w,K1,M1,K2,M2]=iqc_slowtv(v,d,a,D)
% 
%   Adds iqc-terms describing the block
%   w(t)=k(t)*v(t), |k(t)|<= D,  | d(k(t))/dt | <= d.
%   
%   v : allowed to be a vector signal.
%   a : poles of the multiplier in the IQC.   
%
%   See the reference manual for the detail of the IQCs
%   described in this M-file. 
%   
%   defalut :  d = 1, D = 1, a = 1. 
%
% Created cykao@mit.edu, last modified 01/14/2000 

if      nargin<2,
        D=1;
        a=1;
        d=1;
elseif  nargin<3,
        D=1;
        a=1;
elseif  nargin<4,
        D=1;
end

m=size(v,1);
n=length(a);
A=diag(-a);
B1=ones(n,m);
B2=eye(n);
C=eye(n);
D1=zeros(n,m);
D2=zeros(n);

G1=ss(A,B1,C,D1);
G2=ss(A,B2,C,D2);
w=signal(m);
u=signal(n);
y=G1*v;
z=G1*w;
x=G2*u;

v_ext = [y;v];
w_ext = [z+x;w];

K1 = symmetric(size(v_ext,1));
M1 = skew(size(v_ext,1));
IQC1 = v_ext'*D^2*K1*v_ext-w_ext'*K1*w_ext+v_ext'*(M1*w_ext)+(w_ext'*M1')*v_ext;
IQC1 > 0;
K1 > 0;

K2 = symmetric(size(y,1));
if size(y,1)==1,
   IQC2 = y'*d^2*K2*y - u'*K2*u;
   IQC2 > 0;
   K2 > 0;
else
   M2 = skew(size(y,1));
   IQC2 = y'*d^2*K2*y - u'*K2*u + y'*(M2*u) + (u'*M2')*y;
   IQC2 > 0;
   K2 > 0 ;
end


% ---------------------------------------------- % 
%    old version -- obsolete, kept for record    %
% ---------------------------------------------- % 
%
% G1=ss(A,B1,C,D1);
% G2=ss(A,B2,C,D2);
% w=signal(m);
% w2=signal(n);
% y=G1*v;
% z=G1*w;
% x=G2*w2;
% w1=z+x;

% V=[y;v];
% W=[z;w];
% W1=[w1;w];
% X=[x;zeros(size(w,1),1)];
% K=symmetric(size(V,1));
% K>0;
% V'*K*V-W'*K*W>W1'*2*K*X-X'*K*X;
% M=skew(size(V,1));
% V'*M*W+W'*M'*V==-V'*M*X-X'*M'*V;
% w ==iqc_tvscalar(v,D);  
% w1==iqc_tvscalar(y,D);
% w2==iqc_tvscalar(y,d);
