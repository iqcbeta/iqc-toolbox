%
% A test for iqc_polytope_stvp
%


%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu




clear all 
A=[-1.4928 -0.9898 1.1380 -0.3306;
    -0.6451 -0.1645 -0.6841 -0.8436;
    0.8057 0.2895 -2.7960 0.4978;
    0.2316 1.4789 -0.0729 -0.0156];
B=[-0.5465 -0.8542;
    -0.8468 -1.2013;
    -0.2463   -0.1199;
     0.6630   -0.0653];
C=[0.4853 -0.1497 -0.0793 -0.6065;
    -0.5955 -0.4348 1.5352 -1.347];
G=ss(A,B,C,zeros(2,2));
G=G/(norm(G,Inf)-1);
k=1.96;
Delta{1}=k*[1 0;0 1];
Delta{2}=k*[1 0;0 -1];
Delta{3}=k*[-1 0;0 0];
Delta{4}=k*[-1 0;0 -1];
d1=0.1;
d2=0.1;
Omega{1}=[d1 0;0 d2];
Omega{2}=[d1 0;0 -d2];
Omega{3}=[-d1 0;0 d2];
Omega{4}=[-d1 0;0 -d2];
Lambdastruc=[1 0;0 2];
abst_init_iqc;
w=signal(2);
f=signal(2);
v=G*(w+f);
[w1,x,y,z,L]=iqc_polytope_stvp(v,Delta,Omega,Lambdastruc);
w==w1;
gain=iqc_gain_tbx(f,v)
disp('the result was gain=199.66 on Nov 12, 1997');
iqc_bode
iqc_value
La=value_iqc(L)

