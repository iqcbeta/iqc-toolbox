%
% Test of iqc_slope.m and iqc_white.m
%
clear all
G=ss([-0.5 -5;1 0],[1;0],[1 1],0);
Gr=ss(-0.1,9,1,10);
abst_init_iqc;
w=signal;
b=100;
f=iqc_white(1,b,[1,-0.25+2.2*i]);
v=-Gr*G*(f+w);
y=G*(f+w);
a=2.5;
N=2;
alpha=0;
beta=1;
w==iqc_slope(v,a,N,alpha,beta);
gain=iqc_gain_tbx(f,y)
iqc_bode
disp('the result was gain=0.4537 on dec 2, 1997'); 