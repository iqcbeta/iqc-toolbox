function iqc_gain_tbx_test1(n,m,k,cs)
% function iqc_gain_tbx_test1(n,m,k,cs)
%
% testing the function iqc_gain_tbx
% by calculating the L2 induced norm f->y of a
% randomly generated LTI system 
% dx/dt=Ax+Bf, y=Cx+Df
% n=size of x
% m=size of f
% k=size of y
% 
% default n=5,m=2,k=2,cs=1
%
% Written by ameg@mit.edu,  last modified October 13, 1997

if nargin<1, n=5; end   % default arguments
if nargin<2, m=2; end
if nargin<3, k=2; end
if nargin<4,cs=1; end

randn('state',cs);      % set the random generator seed
A=randn(n);             % generate the system
B=randn(n,m);
C=randn(k,n);
D=randn(k,m);

A=A-(1+max(real(eig(A))))*eye(n);  % make the system stable
G=ss(A,B,C,D);

disp('*** Testing iqc software on an H-infinity norm calculation problem')
disp('*')
disp(['*     ||G|| = ' num2str(norm(ss(A,B,C,D),Inf))])
disp(' ')

abst_init_iqc

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

f=signal(m);
iqc_gain_tbx(f,G*f)
iqc_bode
