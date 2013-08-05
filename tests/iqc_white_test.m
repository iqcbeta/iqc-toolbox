%
% iqc_white_test  (Test of iqc_white)
%
clear all
G11=ss(-1,1,1,0);
G12=ss([-0.2 -1;1 0],[1;0],[1 2],0);
G21=ss(-0.5,1,2,0);
G22=ss(-0.1,0.1,1,0);
G=[G11 G12;G21 G22];
b=50;
H=sqrt(b*2/pi);
a=[0.1;0.5;1;-0.1+0.9950i];
abst_init_iqc;
[f,Y,X]=iqc_white(2,b,a);
z=G*H*f;
gain=iqc_gain_tbx(f,z)
display(['answer as b-> inf: ||G||_{H2}=',num2str(norm(G,2))])
iqc_bode

