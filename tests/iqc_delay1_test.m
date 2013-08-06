%
% Test of iqc_delay1
%
clear all
G=tf(0.5,[1 0.25 1]);
abst_init_iqc;

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sdpt3')

w=signal;
f=signal;
v=G*(-w+f);
w==iqc_delay1(v,0.5);
gain=iqc_gain_tbx(f,w)
disp('the result was gain=185.31 on June 22, 1999')
iqc_bode
