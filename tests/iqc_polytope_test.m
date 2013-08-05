%
% A test for iqc_polytope.m (it shows that the toolbox works for the
% case with nonsquare G and Delta
% 

G1=ss([-0.3 -1;1 0],[1;0],[0 0.15],0);
G2=ss(-1,-0.72,1,0);
G=[G1 G2];
Delta{1}=[1;1];
Delta{2}=[1;-1];
Delta{3}=[-1;1];
Delta{4}=[-1;-1];

abst_init_iqc;
w=signal(2);
f=signal(2);
v=G*(w+f);
w==iqc_polytope(v,Delta);
gain=iqc_gain_tbx(f,v)
disp('the result was gain=234.28 on nov 10, 1997');
iqc_bode
