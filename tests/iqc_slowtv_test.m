function gain=iqc_slowtv_test(k1,k2)
% function iqc_tvscalar_test(k1,k2)
% This is a test file of IQC of slow time varying
% varying scalar 
% 
% system: xdot = Ax + B(delta*Cx);
% where : A=[-0.21,-1;1,0]
%         B=[0.8;0]
%         C=[0,1]
%         |delta(t)|<=k2
%         |delta(t)'|<=k1
%
% It has been verified that when k1=1. and k2=0.311, system is stable.
% Also, system is unstable when k2=1 k1>0.46.
% For k2=1, the best k1 this program can get is 0.1545.

if nargin==0
   k1=0.1;
   k2=1;
elseif nargin==1,
   k2=1;
end

if k1<0,
error('parameter k1 should be non-negative')
end

if k2<=0,
error('parameter k2 should be positive')
end

abst_init_iqc;

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

A=[-0.21,-1;1,0];
B=[0.8;0];
C=[0,1];
D=0;
G=ss(A,B,C,D);

w=signal;
f=signal;
v=G*(f+w);
w==iqc_slowtv(v,k1,[1,3,5],k2);

gain=iqc_gain_tbx(f,w);

if ~isempty(gain)
   iqc_bode
end