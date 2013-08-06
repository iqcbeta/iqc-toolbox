%
%test of iqc_dzn_e
%
clear all
s=tf([1 0],1);
P=-1/(s*s+s+1)
k1=3;
k2=0.3;
kphi=k2*1;
G=(s+1-k1/k2)/(s*s+s+1);
abst_init_iqc;

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

w=signal;
f=signal;
v=G*w+f;
w==iqc_dzn_e(v,1.8,4,kphi);
gain=iqc_gain_tbx(f,v)
iqc_bode
display('result was gain=2.32 on dec 5, 1997')
