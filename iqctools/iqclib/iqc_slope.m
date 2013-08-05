function [w,h_0,H,x]=iqc_slope(v,a,N,alpha,beta)
% function [w,h_0,H,x]=iqc_slope(v,a,N,alpha,beta)
%
% Defines the IQC for a slope restricted nonlinearity, i.e.
% 
%       (beta*v-w)'*[(h_0-H)*(w-alpha*v)]>0,
%
% where 
%       N
% H(s)=sum x_k /(s+a)^{k+1}
%      k=0
%
% Here a is any nonzero real value and N>=0.
%
% The following constraints need to be satisfied
%
%                                         /
%                                        | 0, t>0 if a>0
%  1. h(t)=sum_{k=0}^N x_k t^ke^{-at}/k!>|
%                                        | 0, t<0 if a<0
%                                         \
%      N
%  2. sum x_k/a^{k+1}
%     k=0
%
% Default values: a=1
%                 N=0
%                 alpha=0
%                 beta=1

%
% Work by Ulf Jonsson in Oct and Nov of 1997
%

if nargin<1
   error('input must be specified')
end
if size(v,1)~=1, error('only scalar signals are allowed'), end
if nargin<2, a=1; end
if nargin<3, N=0; end
if nargin<4, alpha=0; end
if nargin<5, beta=1; end

sn=sign(a);
a=abs(a);
s=tf([1 0],1);
h_0=symmetric;
h_0>0;
if N==0
  x=symmetric;
  x>0;
  x<h_0*a;
  H=x*abst(1/(s+a));
else
  x=rectangular(1,N+1);
  A=-a*eye(N+1);
  E=eye(N);
  A(2:N+1,1:N)=A(2:N+1,1:N)+E;
  B=[1;zeros(N,1)];
  C=eye(N+1);
  D=zeros(N+1,1);
  H=x*ss(A,B,C,D);
  x*(-A\B)<h_0;
  Vbinom=[];
  FM=[1];
  for k=1:N
    Vbinom=[gamma(N+1)/(gamma(k+1)*gamma(N-k+1)) Vbinom];
    FM=[FM;(-1)^k/gamma(k+1)];
  end
  FM=diag(FM);
  as=Vbinom; 
  A=[zeros(N-1,1) eye(N-1);-as];
  B=[zeros(N-1,1);1];
  C=[eye(N);-as];
  D=[0;B];
  as=[as';1];
  M=[];
  for k=0:N
    Mk=zeros(N+1,1);
    for l=0:min(2*k,N)
      if 2*k-l<=N
        Mk(l+1)=as(2*k-l+1)*(-1)^l;
      end
    end
    M=[M Mk];
  end
  F=x*(FM/M)*ss(A,B,C,D);
  F>0;
end
w=signal;
if sn>0
   (w-alpha*v)'*(h_0-H)*(beta*v-w)>0;
else
   ((h_0-H)*(w-alpha*v))'*(beta*v-w)>0;
end




