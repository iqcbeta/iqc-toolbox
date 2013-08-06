%
% A test program for iqc_popov_vect.m
% 
clear all
A=[-0.37 0.20 0.15;...
   -0.24 -0.65 0.51;...
    0.09 -0.53 -0.60];
B=[-0.14 0;...
   0.11 -0.10;...
    0   -0.83];
C=[0.15   0    0;...
   0     0.8   0.4];
D=zeros(2,2);
G=ss(A,B,C,D);
s=tf([1 0],1);
abst_init_iqc;

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')


w=signal(2);
f=signal(2);
v=G*w+1/(s+1)*f;
w==iqc_popov_vect(v);
w==iqc_tvscalar(v,1);
gain=iqc_gain_tbx(f,v)
