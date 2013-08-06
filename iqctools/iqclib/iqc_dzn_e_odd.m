function [w,h_0,H,F,xp,xm]=iqc_dzn_e_odd(v,a,N,k_dzn)
% function [w,h_0,H,F,xp,xm]=iqc_dzn_e_odd(v,a,N,k_dzn)
%
% Defines IQCs for a integrator encapsulated in an odd deadzone
% nonlinearity with slope in [0,k_dzn].
% The operator is defined by
%
%       dz/dt=dzn(v-z),  z(0)
%           w=dzn(v-z)
%
% where dzn denotes the deadzone nonlinearity.
%
% The IQCs are defined as
%
%       w'*(h_0-H)*(v-w/k_dzn)+w'*Fw>0
%
% where F=(H(s)-H(0))/s, and
%       N
% H(s)=sum x_k k!/(s+a)^{k+1}
%      k=0
%
% Here a is any nonzero real value and N>=0.
%
% The following constraint need to be satisfied
%
%  ||h||_1<h_0, where h is the weighting function corresponding to H
%
% Default values: a=1
%                 N=0
%                 k_dzn=1
%
%
% Work by Ulf Jonsson in December of 1997
% Last modified by cmj on 2013/5/2

if nargin<1
    disp_str(1)
end
if size(v,1)~=1
    disp_str(8,'scalar')
end
if nargin<2
    a=1;
end
if nargin<3
    N=0;
end
if nargin<4
    k_dzn=1;
end

global ABST

switch ABST.systemtype
    case 'continuous'
        sn=sign(a);
        a=abs(a);
        s=tf([1 0],1);
        h_0=symmetric;
        h_0>0; %#ok<*VUNUS>
        if N==0
            xp=symmetric;
            xm=symmetric;
            xp>0;
            xm>0;
            xp+xm<h_0*a;
            H=(xp-xm)*ss(-a,1,1,0);
            F=(xp-xm)*ss(-a,-1/a,1,0);
        else
            xp=rectangular(1,N+1);
            xm=rectangular(1,N+1);
            A=-a*eye(N+1);
            E=eye(N);
            A(2:N+1,1:N)=A(2:N+1,1:N)+E;
            B=[1;zeros(N,1)];
            C=eye(N+1);
            D=zeros(N+1,1);
            H=(xp-xm)*ss(A,B,C,D);
            F=(xp-xm)*ss(A,A\B,C,D);
            (xp+xm)*(-A\B)<h_0;
            Vbinom=[];
            FM=1;
            for k=1:N
                Vbinom=[gamma(N+1)/(gamma(k+1)*gamma(N-k+1)) Vbinom]; %#ok<*AGROW>
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
            Fp=xp*(FM/M)*ss(A,B,C,D);
            Fp>0;
            Fm=xm*(FM/M)*ss(A,B,C,D);
            Fm>0;
        end
        w=signal;
        if sn>0
            w'*(h_0-H)*(v-(1/k_dzn)*w)+w'*F*w>0;
        else
            ((h_0-H)*w)'*(v-(1/k_dzn)*w)+(F*w)'*w>0;
        end
    case 'discrete'
        disp_str(70,'iqc_dzn_e_odd','discrete')
end
