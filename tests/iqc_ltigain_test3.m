clc
a=-0.3;
b=0.8;
k=2.5;
Td=0.86;
A=[-(a+b*k*Td), -b*k; 1 0];
B=[-2, -2*b; 0, 0];
C=[a+b*k*Td, b*k; k*Td, k];
D=[1, 2*b;0, 1];
G=ss(A,B,C,D);
abst_init_iqc;
w=signal(2);
f=signal(2);
v=G*(f+w);
% [w1,p,q,x,y,z]=iqc_ltigain(v);
% w==0.22*w1;
w==0.22*iqc_ltigain(v);
gain=iqc_gain_tbx(f,v)
disp('the result was gain=11.1456 on May 29, 2013');