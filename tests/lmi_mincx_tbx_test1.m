function lmi_mincx_tbx_test1(n,m,k,cs)
% function lmi_mincx_tbx_test1(n,m,k,cs)
%
% testing the function lmi_mincx_tbx
% by calculating the square if the L2 induced norm f->y of a
% randomly generated LTI system
% dx/dt=Ax+Bu, y=Cx+Df
% n=size of x
% m=size of f
% k=size of y
% 
% solves the system of LMI's 
%   [P*A+A'*P+C'*C    P*B+C'*D;
%     B'*P+D'*C       D'*D-y*I(m) ] <0,    y->min
%
% default n=10,m=2,k=2,cs=1
%
% Written by ameg@mit.edu,  last modified October 13, 1997

%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu



if nargin<1, n=10; end   % default arguments
if nargin<2, m=2; end
if nargin<3, k=2; end
if nargin<4,cs=1; end

randn('state',cs);      % set the random generator seed
A=randn(n);             % generate the system
B=randn(n,m);
C=randn(k,n);
D=randn(k,m);

A=A-(1+max(real(eig(A))))*eye(n);  % make the system stable

disp('*** Testing lmi software on an H-infinity norm calculation problem')
disp('*')
disp(['*     ||G||^2 = ' num2str(norm(ss(A,B,C,D),Inf)^2)])
disp(' ')

abst_init_lmi;
p=symmetric(n);
y=symmetric;
a=abst(A);
b=abst(B);
c=abst(C);
d=abst(D);
[p*a+a'*p+c'*c p*b+c'*d;b'*p+d'*c d'*d-y*II(m)]<0;
lmi_mincx_tbx(y);
Y=value_iqc(y)

