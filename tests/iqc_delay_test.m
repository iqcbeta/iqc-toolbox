% 
% A test program for function iqc_delay.m
%
clear all
G=tf(0.5,[1 0.25 1]);
abst_init_iqc;
w=signal;
f=signal;
v=G*(-w+f);
w==iqc_delay(v,0.5,[1,2]);
gain=iqc_gain_tbx(f,w)
disp('the result was gain=1.87 on June 22, 1999')
iqc_bode
