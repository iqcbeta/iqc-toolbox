%
% Another test program for iqc_polytope_stvp
% 
clear all
A=[-3 0 0 0;
   0 -0.01 1 0;
   0 -1 -0.01 0;
   0  0   0  -1];
B=0.002*[0 4;-1 3;4 5;8 14];
C=0.1*[-3 15 5 -5;
       15 20 5 -5];
G=ss(A,B,C,zeros(2,2));
Delta{1}=[1 0;0 1];
Delta{2}=[1 0;0 -1];
Delta{3}=[-1 0;0 0];
Delta{4}=[-1 0;0 -1];
Omega{1}=[0.1 0;0 0.1];
Omega{2}=[0.1 0;0 -0.1];
Omega{3}=[-0.1 0;0 0.1];
Omega{4}=[-0.1 0;0 -0.1];
Lambdastruc=[1 0;0 2];
abst_init_iqc;

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sdpt3')

w=signal(2);
f=signal(2);
v=G*(w+f);
w==iqc_polytope_stvp(v,Delta,Omega,Lambdastruc);
gain=iqc_gain_tbx(f,v)