disp('***********************************************************')
disp('*       This is a global test of the iqc_beta package     *')
disp('*  Written by ameg@mit.edu, last modified October 13,1997 *')
disp('***********************************************************')


%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu




disp(' ')
disp(' ')
abst_init_iqc
lmitbx_options([0 0 0 0 1])


disp(' ')
disp(' ')
disp('*** iqc_gain_tbx_test1 ...')
n=5;m=2;k=2;
randn_state=randn('state');
A=randn(n);             % generate the system
B=randn(n,m);
C=randn(k,n);
D=randn(k,m);
A=A-(1+max(real(eig(A))))*eye(n);  % make the system stable
G=ss(A,B,C,D);
abst_init_iqc
f=signal(m);
g=iqc_gain_tbx(f,G*f);
ga=norm(G,Inf);
disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if (g<ga)|(g>1.05*ga),
   error('iqc_beta is corrupted')
end
disp(' ')
disp('*** iqc_value_test1 ...')
s=tf([1 0],1);
G=(s-1)/((s+2)*(s+2));
H=(s-2)/((s+1)*(s+1));
ga=2;
abst_init_iqc;
lmitbx_options([0 0 0 0 1]);
w=signal(2);
v=[G 1]*w;
x=symmetric;
x>0;
sigma=[v' 0]*x*[0 H;0 0]*[0;w(1)]-w(1)'*x*w(1);
sigma>0;
g=iqc_gain_tbx(w(2),v);
iqc_value;
[X,Sigma,W,V]=value_iqc(x,sigma,w,v);
d=[V' [0;0]]*X*[0 H;0 0]*[[0 0];W(1,:)]-W(1,:)'*X*W(1,:)-Sigma;
er=norm(d,Inf);
disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if (g<ga)|(g>1.05*ga)|(er>1e-10),
   error('iqc_beta is corrupted')
end

disp(' ')
disp('*** iqc_get_mlmi_test ...')
s=tf([1 0],1);
G=(s-1)/((s+2)*(s+2));
abst_init_iqc;
lmitbx_options([0 0 0 0 1]);
w=signal;
v=G*w;
gain=iqc_gain_tbx(w,v);
[P,A,B,Sigma,MainLMI]=iqc_get_mlmi;
[AG,BG,CG,DG]=ssdata(G);
T=[CG,DG;zeros(size(DG,2),size(CG,2)),eye(size(DG,2))];
Sigma_verify=T'*[1,0;0,-gain*gain]*T;
if max(svd(Sigma_verify-Sigma))>1e-3
   error('iqc_beta is corrupted')
end

disp(' ')
disp('*** iqc_ltigain_test1 ...')
randn_state=rand('state');
g=rand(1,4);                   % coefficients of g(s)
h=rand(1,4);                   % coefficients of h(s)
sg=tf(g(1:2),[1 g(3:4)]);      % CS Toolbox representation of g
sh=tf(h(1:2),[1 h(3:4)]);      % CS Toolbox representation of h
ga=sqrt(1+max(norm(sg+sh,inf),norm(sg-sh,inf))^2);  % analytical gamma
G=[0 sg;0 0];
H=[sh;1];
abst_init_iqc
w=signal(2);
f=signal;
v=G*w+H*f;
w==iqc_ltigain(v);
g=iqc_gain_tbx(f,v);
disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if (g<ga)|(g>1.05*ga),
   error('iqc_beta is corrupted')
end

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
 g=iqc_gain_tbx(v,y);
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
g=iqc_gain_tbx(f,v);
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
g=iqc_gain_tbx(f,w);
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
g=iqc_gain_tbx(f,w);
disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if (abs(g-ga)>ga/20)|(g<ga),
   error('iqc_beta is corrupted')
end

disp(' ')
disp('*** iqc_white_test ...')
G11=ss(-1,1,1,0);
G12=ss([-0.2 -1;1 0],[1;0],[1 2],0);
G21=ss(-0.5,1,2,0);
G22=ss(-0.1,0.1,1,0);
G=[G11 G12;G21 G22];
b=50;
H=sqrt(b*2/pi);
a=[0.1;0.5;1;-0.1+0.9950i];
ga=norm(G,2);
abst_init_iqc;
[f,Y,X]=iqc_white(2,b,a);
z=G*H*f;
g=iqc_gain_tbx(f,z);
disp(['   Upper bound: ' sft(g,10)]);
%if (abs(g-ga)>ga/20),
%   error('iqc_beta is corrupted')
%end

disp(' ')
disp('*** lmi_mincx_tbx_test1 ...')
n=5;m=2;k=2;
randn_state=randn('state');
A=randn(n);             % generate the system
B=randn(n,m);
C=randn(k,n);
D=randn(k,m);
A=A-(1+max(real(eig(A))))*eye(n);  % make the system stable
ga=norm(ss(A,B,C,D),Inf);
abst_init_lmi;
p=symmetric(n);
y=symmetric;
a=abst(A);
b=abst(B);
c=abst(C);
d=abst(D);
[p*a+a'*p+c'*c p*b+c'*d;b'*p+d'*c d'*d-y*II(m)]<0;
g=sqrt(lmi_mincx_tbx(y));
[Y,P]=value_iqc(y,p);
lm=min(eig(-[P*A+A'*P+C'*C P*B+C'*D;B'*P+D'*C D'*D-Y*eye(m)]));
disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if (abs(g-ga)>ga/20)|(lm<0),
   error('iqc_beta is corrupted')
end

disp(' ')
disp('*** lmi_mincx_tbx_test2 ...')
n=5; 
randn_state=randn('state');
A=randn(n);             % generate the data
B=randn(n);
d=eig(0.5*(A'+A-B'-B));
ga=sum(d(find(d>0)))+trace(B);
abst_init_lmi;
x=symmetric(n);
x>A;
x>B;
g=lmi_mincx_tbx(trace(x));
X=value_iqc(x);
lm=min([eig(X-0.5*(A+A'));eig(X-0.5*(B+B'))]);
disp(['      Lower bound: ' sft(ga,10) ';  Upper bound: ' sft(g,10)]);
if (abs(g-ga)>ga/20)|(lm<0)|(abs(g-trace(X))>0.01),
   error('iqc_beta is corrupted')
end

disp(' ')
disp('*** fdlmi_mincx_tbx_test ...')
s=tf([1,0],1);
abst_init_fdlmi;
lmitbx_options([0 0 0 0 1]);
G0=1/(s*s+s+1);
G1=ss(-1,1,1,0);
G2=ss(-0.7,1,1,0);
x0=symmetric;
x1=symmetric;
x2=symmetric;
y =symmetric;
y>0;
Ga=x0+x1*G1+x2*G2;
H = (G0-Ga);
[y,H;H',y]>0;
fdlmi_mincx_tbx(y);
if abs(value_iqc(y)-0.3012)>0.01,
   error('iqc_beta is corrupted')
end


disp(' ')
disp('*** iqc_cdelay_test ...')
G=ss([-1 -2;1 0],[1;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=G*(w+f);
w==iqc_cdelay(v,0.5,[1 2]);
gain=iqc_gain_tbx(f,v);
if abs(gain-9.5)>0.2,
   error('iqc_beta is corrupted')
end

disp('*** iqc_ltiunmod_test ...')
s=tf([1 0],1);
M=1/1.5*[s/(s+100)    0;
         s/(s+80)     0;
        10/(s+10)   8/(s+8)];
abst_init_iqc
f=signal(2);
w1=signal;
w2=signal;
a=4;
v=M*(f+[w1;w2]);
[waux1,x]=iqc_ltiunmod(v(1),a,1);
[waux2,x]=iqc_ltiunmod(v(2:3),a,1);
w1==waux1;
w2==waux2;
g=iqc_gain_tbx(f,v);
if abs(g-2.8294)>0.05*2.8294,
   error('iqc_beta is corrupted')
end


disp('*** iqc_sector_test ...')
G=ss([-0.2 -1;1 0],[1;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=-G*(w+f);
w==iqc_sector(v,0,1);
gain=iqc_gain_tbx(f,v);
if abs(gain-25.97)>0.5,
   error('iqc_beta is corrupted')
end

disp('*** iqc_popov_vect_test ...')
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
w=signal(2);
f=signal(2);
v=G*w+1/(s+1)*f;
w==iqc_popov_vect(v);
w==iqc_tvscalar(v,1);
gain=iqc_gain_tbx(f,v);
if abs(gain-7.7680)>0.3,
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
gain=iqc_gain_tbx(f,v);
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
gain=iqc_gain_tbx(f,v);
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
gain=iqc_gain_tbx(f,v);
if abs(gain-199.66)>5,
   error('iqc_beta is corrupted')
end

disp('*** iqc_window_test ...')
G=ss([-0.3 -100;1 0],[1.4;0],[1 1],0);
abst_init_iqc;
w=signal;
f=signal;
v=G*(w+f);
w==iqc_window(v,0.5);
gain=iqc_gain_tbx(f,v);
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
gain=iqc_gain_tbx(f,y);
if abs(gain-29.4806)>1,
   error('iqc_beta is corrupted')
end


disp('*** iqc_delay_test ...')
G=tf(0.5,[1 0.25 1]);
abst_init_iqc;
w=signal;
f=signal;
v=G*(-w+f);
w==iqc_delay(v,0.5,[1,2]);
gain=iqc_gain_tbx(f,w);
if abs(gain-1.87)>0.5,
   error('iqc_beta is corrupted')
end


disp('*** iqc_delay1_test ...')
G=tf(0.5,[1 0.25 1]);
abst_init_iqc;
w=signal;
f=signal;
v=G*(-w+f);
w==iqc_delay1(v,0.5);
gain=iqc_gain_tbx(f,w);
if abs(gain-185.31)>2,
   error('iqc_beta is corrupted')
end


disp('*** iqc_slowtv_test ...')
abst_init_iqc;
A=[-0.21,-1;1,0];
B=[0.8;0];
C=[0,1];
D=0;
G=ss(A,B,C,D);
w=signal;
f=signal;
v=G*(f+w);
w==iqc_slowtv(v,0.1,[1,3,5],1);
gain=iqc_gain_tbx(f,w);
if abs(gain-25.7744)>1,
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
gain=iqc_gain_tbx(f,w);
if abs(gain-22.6304)>1,
   error('iqc_beta is corrupted')
end



% disp('*** GUI test1 ...')
% g=iqc_gui('test1');
% ga=1;
% if isempty(g)|(abs(g-ga)>ga/1000),
%    error('iqc_beta is corrupted')
% end
% 
% disp('*** GUI test2 ...')
% g=iqc_gui('test2');
% ga=1.8474;
% if isempty(g)|(abs(g-ga)>ga/1000),
%    error('iqc_beta is corrupted')
% end
% 
% disp('*** GUI test3 ...')
% g=iqc_gui('test3');
% ga=0.6903;
% if isempty(g)|(abs(g-ga)>ga/1000),
%    error('iqc_beta is corrupted')
% end
% 
% disp('*** GUI test4 ...')
% g=iqc_gui('test4');
% ga=0.1429;
% if isempty(g)|(abs(g-ga)>ga/1000),
%    error('iqc_beta is corrupted')
% end

disp(' ')
disp(' ')
disp('***********TESTING OK!!!')
