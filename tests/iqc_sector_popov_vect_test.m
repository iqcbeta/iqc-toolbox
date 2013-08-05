%
% A test of iqc_sector and iqc_popov_vect
%
% First iqc_sectorNL and then iqc_sectorNL+iqc_mvpopov
clear all
G=ss([-0.2 -1;1 0],[1;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=-G*(w+f);
w==iqc_sector(v,0,1);
gain=iqc_gain_tbx(f,v)
disp('The answer was gain=25.97 on nov 10, 1997')

G=ss([-0.2 -1;1 0],[1;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=-G*(w+f);
w==iqc_sector(v,0,1);
w==iqc_popov_vect(v,'0');
gain=iqc_gain_tbx(f,v)
disp('The answer was gain=7.08 on nov 10, 1997')
iqc_bode
