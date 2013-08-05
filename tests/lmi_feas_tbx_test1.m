function lmi_feas_tbx_test1(n,cs)
% function lmi_feas_tbx_test1(n,cs)
%
% testing the function lmi_feas_tbx
% on the problem X+Y=A, lmax(X,Y)->min 
% (A is symmetric of size nxn)
%
% Written by ameg@mit.edu,  last modified October 13, 1997
if nargin<1, n=5; end   % default arguments
if nargin<2, cs=1; end

randn('state',cs);      % set the random generator seed
A=randn(n);             % generate the system
A=A+A';

disp('*** Testing lmi software on an eigenvalue calculation problem')
disp('*')
disp(['*     answer should be = ' num2str(0.5*max(eig(A)))])
disp(' ')

abst_init_lmi;
x=symmetric(n);
x<0;
A<x;
lmi_mincx_tbx

