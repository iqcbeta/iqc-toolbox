G=-ss([-0.1 -1;1 0],[1;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=G*(w+f);
w==iqc_popov(v,0,1);
gain=iqc_gain_tbx(f,v)
disp('the result was gain=14.1492 on May 29, 2013');