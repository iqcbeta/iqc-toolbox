%
% Test of iqc_intdelaydiff.m
%

G=ss([-0.3 -100;1 0],[1.4;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=G*(w+f);
w==iqc_window(v,0.5);
gain=iqc_gain_tbx(f,v)
display('result was gain=627.60 on nov 10, 1997')
iqc_bode