function iqc_ltvnorm_test1(n,m,k,cs)
% function iqc_ltvnorm_test1(n,m,k,cs)
%
% testing "iqc_ltvnorm": iqc description of
% a norm-bounded LTV uncertainty
% 
% applies "iqc_norm" to the system v=Gw+f,
% where v is n-dimensional,
% w is m-dimensional, and a stable transfer matrix G=G(s)
% is generated randomly
%
% the answer obtained using iqc_gain_tbx is
% compared with the analytical answer, which states that
% the square of the worst case induced norm is
% k(k+1)/[k-(k+1)r^2], where r is the H-infinity norm of G
% and k=r/[1-r].
%
% Written by ameg@mit.edu,  last modified October 13, 1997

if nargin<1, n=5; end
if nargin<2, m=3; end
if nargin<3, k=2; end
if nargin<4, cs=1; end

randn('state',cs);
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

k1=hinf/(1-hinf);      % calculating the true gamma
gain=sqrt(k1*(k1+1)/(k1-(k1+1)*(hinf^2)));

disp('******* testing iqc_gain_tbx ***********')
disp('****************************************')
disp(['ANALYTICAL ANSWER = ' num2str(gain)])

abst_init_iqc

setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sdpt3')

w=signal(m);
f=signal(k);
v=G*w+f;
w==iqc_ltvnorm(v,m);
gain=iqc_gain_tbx(f,v)
iqc_bode