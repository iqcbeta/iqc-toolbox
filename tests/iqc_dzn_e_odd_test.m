%
% test of iqc_dzn_e_odd
%
% The example from the encapsulation paper

clear all
s=tf([1 0],1);
P=-1/(s*s+s+1)
k1=3;
k2=0.3;
kphi=k2*1;
Go=(s+1-k1/k2)/(s*s+s+1);
G22=Go;
a=1;
G21=(k1*s+k2)/(kphi*(s+a))
G12=(s+a)*P;
abst_init_iqc;

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sdpt3')

w1=signal;
w2=signal;
f=signal;
v1=G12*w2;
v2=G22*w2+G21*w1+f;
T0=0.27
[w11,X,x]=iqc_window(v1,T0,[]);
w1==w11;
w2==iqc_dzn_e_odd(v2,1.8,0,kphi);
gain=iqc_gain_tbx(f,v2)
iqc_bode
display('result was gain=98.69 on dec 5, 1997')
