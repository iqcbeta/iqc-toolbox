%
% Test of iqc_cdelay.m 
%
% clear all
clc
G=ss([-1 -2;1 0],[1;0],[1 1],0);

abst_init_iqc;

setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sdpt3')

w=signal;
f=signal;
v=G*(w+f);
w==iqc_cdelay(v,0.5,[1 2]);
gain=iqc_gain_tbx(f,v)
disp('the result was gain=9.50 on nov 10, 1997');
iqc_bode


