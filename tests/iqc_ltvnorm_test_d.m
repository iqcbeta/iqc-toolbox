function iqc_ltvnorm_test_d(n,m,k)

if nargin<1, n=5; end
if nargin<2, m=3; end
if nargin<3, k=2; end

A=randn(n);
A=A-(0.1+max(real(eig(A))))*eye(n);
B=randn(n,m);
C=randn(k,n);
D=randn(k,m);

G=ss(A,B,C,D);         % G(s) is defined now as a CS Toolbox LTI

G=c2d(G,1);

[hinf,wo]=norm(G,inf); % H-infinity norm of G using CS Toolbox
kk=hinf+3;             % tuning coefficient
G=G/kk;                % making H-infinity norm of G less than 1
hinf=hinf/kk;          % correcting the H-infinity norm

k1=hinf/(1-hinf);      % calculating the true gamma
gain=sqrt(k1*(k1+1)/(k1-(k1+1)*(hinf^2)));

disp('******* testing iqc_gain_tbx ***********')
disp('****************************************')
disp(['ANALYTICAL ANSWER = ' num2str(gain)])

abst_init_iqc(1)

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sdpt3')

w=signal(m);
f=signal(k);
v=G*w+f;
w==iqc_ltvnorm(v,m);
gain=iqc_gain_tbx(f,v)
iqc_bode(0,pi)