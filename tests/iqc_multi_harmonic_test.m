%
% iqc_multi_harmonic_test
%
clear all
n=2;
I=eye(n);
O=zeros(n,n);
A0=[-3 0;0 -6];
A1=[0 1;2 0];
B1=[1 0;0 -1];
A2=[0 2;1 0];
B2=[-1 1;1 1];
G=ss(A0,[A1 B1 A2 B2 I],I,[O O O O O]);
abst_init_iqc;
w1=signal(4*n);
w2=signal(n);
f=signal(n);
v=G*[w1;w2]+f; 
w1==iqc_multi_harmonic(v,[1 2],2);
w2==iqc_ltvnorm(v,n,0.4);
gain=iqc_gain_tbx(f,v)

disp('the result was gain=74.9938 on June 22, 1999');

