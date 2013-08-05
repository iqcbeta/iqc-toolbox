%
% iqc_domharmonic_test.m
%
clear all
s=tf([1 0],1);
G=100/(s*s+0.1*s+100);
abst_init_iqc;
%%[f,M,x,d]=iqc_domharmonic(1,5,[],9,2);
[f,M,x,d]=iqc_domharmonic(1,5,[],9,1);
z=G*f;
gain=iqc_gain_tbx(f,z)
disp('the result was gain=70.95 on June 22, 1999');
iqc_bode
