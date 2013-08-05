%
% Test of odd slope restricted nonlinearity
%

clear all
s=tf([1 0],1);
G=2*20*10*(s*s-0.5*s+1)/((s+10)*(s+20)*(s+30));

abst_init_iqc;
w=signal;
f=signal;
v=G*(-w+f);
a=-1;
N=1;
alpha=0;
beta=1;
[w1,h_0,H,xp,xm]=iqc_slope_odd(v,a,N,alpha,beta);
w==w1;
gain=iqc_gain_tbx(f,v)
iqc_bode;
disp('the result was gain=28.99 on nov 10, 1997');
