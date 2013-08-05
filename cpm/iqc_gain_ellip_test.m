abst_init_iqc


disp(' ')
disp('*** iqc_ltigain_test2 ...')
rand_state=rand('state');
 cg=rand(1,4);                   % coefficients of g(s)
 ch=rand(1,4);                   % coefficients of h(s)
 g=tf(cg(1:2),[1 cg(3:4)]);      % CS Toolbox representation of g
 h=tf(ch(1:2),[1 ch(3:4)]);      % CS Toolbox representation of h
 ga=max(norm(g+h,inf),norm(g-h,inf));  % analytical gamma
 abst_init_iqc
 w=signal;
 v=signal;
 y=g*w+h*v;
 w==iqc_ltigain(v);
 g=iqc_gain_ellip(v,y);
 disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if abs(g-ga)>ga/20,
   error('iqc_beta is corrupted')
end



disp(' ')
disp('*** iqc_ltvnorm_test1 ...')
n=5;m=3;k=2;
randn_state=randn('state');
A=randn(n);
A=A-(0.1+max(real(eig(A))))*eye(n);
B=randn(n,m);
C=randn(k,n);
D=randn(k,m);
G=ss(A,B,C,D);         % G(s) is defined now as a CS Toolbox LTI
[hinf,wo]=norm(G,inf); % H-infinity norm of G using CS Toolbox
kk=hinf+3;             % tuning coefficient
G=G/kk;                % making H-infinity norm of G less than 1
hinf=hinf/kk;          % correcting the H-infinity norm
k1=hinf/(1-hinf);          % calculating the true gamma
ga=sqrt(k1*(k1+1)/(k1-(k1+1)*(hinf^2)));
abst_init_iqc
w=signal(m);
f=signal(k);
v=G*w+f;
w==iqc_ltvnorm(v,m);
g=iqc_gain_ellip(f,v);
disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if abs(g-ga)>ga/20,
   error('iqc_beta is corrupted')
end



disp(' ')
disp('*** iqc_popov_test1 ...')
G=tf([-1 1],[1 2 5]);
G=G/(norm(G,inf)+1);
[m,p,ww]=bode(G);
ww=linspace(0.001,max(ww)/3,1000);
g=freqresp(G,ww);
g=squeeze(g);
lbw=abs(g./(1+g)).*(real(g)>-1)+abs(g./imag(g)).*(real(g)<=-1);
ga=max(lbw);
abst_init_iqc
f=signal;
w=signal;
v=G*(f-w);
w==iqc_popov(v,[0 1]);
g=iqc_gain_ellip(f,w);
disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if abs(g-ga)>ga/20,
   error('iqc_beta is corrupted')
end



disp(' ')
disp('*** iqc_monotonic_test1 ...')
G=tf([-1 1],[1 2 5]);
G=G/(norm(G,inf)+1);
a=1;
[m,p,ww]=bode(G);
ww=linspace(0.001,max(ww)/3,1000);
g=freqresp(G,ww);
g=squeeze(g);
lbw=abs(g./(1+g)).*(real(g)>-1)+abs(g./imag(g)).*(real(g)<=-1);
ga=max(lbw);
abst_init_iqc
f=signal;
w=signal;
v=G*(f-w);
w==iqc_monotonic(v,a);
g=iqc_gain_ellip(f,w);
disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if (abs(g-ga)>ga/20)|(g<ga),
   error('iqc_beta is corrupted')
end


disp('*** iqc_sector_test ...')
G=ss([-0.2 -1;1 0],[1;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=-G*(w+f);
w==iqc_sector(v,0,1);
gain=iqc_gain_ellip(f,v);
if abs(gain-25.97)>0.5,
   error('iqc_beta is corrupted')
end



disp('*** iqc_sector_popov_vect_test ...')
G=ss([-0.1 -1;1 0],[1;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=-G*(w+f);
w==iqc_sector(v,0,1);
w==iqc_popov_vect(v,'0');
gain=iqc_gain_ellip(f,v);
if abs(gain-14.15)>0.7,
   error('iqc_beta is corrupted')
end


disp('*** iqc_polytope_test ...')
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
gain=iqc_gain_ellip(f,v);
if abs(gain-234.28)>5,
   error('iqc_beta is corrupted')
end



disp('*** iqc_polytope_stvp_test ...')
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
w==iqc_polytope_stvp(v,Delta,Omega,Lambdastruc);
gain=iqc_gain_ellip(f,v);
if abs(gain-199.66)>5,
   error('iqc_beta is corrupted')
end


disp('*** iqc_sector_popov_vect_test ...')
G=ss([-0.1 -1;1 0],[1;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=-G*(w+f);
w==iqc_sector(v,0,1);
w==iqc_popov_vect(v,'0');
gain=iqc_gain_ellip(f,v);
if abs(gain-14.15)>0.7,
   error('iqc_beta is corrupted')
end


disp('*** iqc_window_test ...')
G=ss([-0.3 -100;1 0],[1.4;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=G*(w+f);
w==iqc_window(v,0.5);
gain=iqc_gain_ellip(f,v);
if abs(gain-627.60)>10,
   error('iqc_beta is corrupted')
end


disp('*** iqc_harmonic_test ...')
a=1;
b=2.5;
w0=3;
s=tf([1 0],1);
abst_init_iqc;
f=signal;
w=signal;   % w=cos(w0*t)*x
y=(1/(s*s+s+a))*(f-b*w);
w==iqc_harmonic(y,w0,[0.7 0.8 5]);
gain=iqc_gain_ellip(f,y);
if abs(gain-29.4806)>1,
   error('iqc_beta is corrupted')
end


disp('*** iqc_delay1_test ...')
G=tf(0.5,[1 0.25 1]);
abst_init_iqc;
w=signal;
f=signal;
v=G*(-w+f);
w==iqc_delay1(v,0.5);
gain=iqc_gain_ellip(f,w);
if abs(gain-185.31)>2,
   error('iqc_beta is corrupted')
end


disp('*** iqc_tvscalar_test ...')
abst_init_iqc;
A=[-0.21,-1;1,0];
B=[0.8;0];
C=[0,1];
D=0;
G=ss(A,B,C,D);
w=signal;
f=signal;
v=G*(f+w);
w==iqc_tvscalar(v,0.25);
gain=iqc_gain_ellip(f,w);
if abs(gain-22.6304)>1,
   error('iqc_beta is corrupted')
end
