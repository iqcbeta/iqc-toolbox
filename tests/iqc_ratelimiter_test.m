A=[-5 1;2 -3];
B=[0;1];
C=[0 1];
D=0;
G=ss(A,B,C,D);

abst_init_iqc

f=signal;
w=signal;
z=G*(f+w);
w==iqc_ratelimiter(z,[1 1],1);

gain=iqc_gain_tbx(f,z)
disp('the result was gain=0.6903 on May 29, 2013');