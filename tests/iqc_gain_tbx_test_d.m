function iqc_gain_tbx_test_d(n,m,k)

if nargin<1
    n=5;
end   % default arguments
if nargin<2
    m=2;
end
if nargin<3
    k=2;
end

A=randn(n);             % generate the system
B=randn(n,m);
C=randn(k,n);
D=randn(k,m);

A=A-(1+max(real(eig(A))))*eye(n);  % make the system stable
G=ss(A,B,C,D);
G=c2d(G,1);
A=G.a;
B=G.b;
C=G.c;
D=G.d;

disp('*** Testing iqc software on an H-infinity norm calculation problem')
disp('*')
disp(['*     ||G|| = ' num2str(norm(ss(A,B,C,D,-1),Inf))])
disp(' ')

abst_init_iqc(1)

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

f=signal(m);
iqc_gain_tbx(f,G*f)
iqc_bode(0,pi)
