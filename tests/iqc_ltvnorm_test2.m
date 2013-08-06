clc
s=tf([1,0],1);
G=0.8/(s*s+0.21*s+1);
abst_init_iqc
f=signal;
w=signal;
v=G*(f+w);
w==iqc_ltvnorm(v,1,0.26);
gain=iqc_gain_tbx(f,v)
disp('the result was gain=954.2745 on May 29, 2013');
