function gain=iqc_tvscalar_test(k)
% function iqc_tvscalar_test(k)
% This is a test file of IQC of arbitrary fast time
% varying scalar
%
% system: xdot = Ax + B(d*Cx);
% where : A=[-0.21,-1;1,0]
%         B=[0.8;0]
%         C=[0,1]
%         |d(t)|<=k,
%
% It has been verified that system is stable when k<=0.26.

if nargin==0
   k=0.26;
end

if k<0,
error('parameter d should be non-negative')
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
w==iqc_tvscalar(v,k);
gain=iqc_gain_tbx(f,w);
if ~isempty(gain)
   iqc_bode
end