function lmi_feas_tbx_test2(n,cs)
% function lmi_feas_tbx_test2(n,cs)
%
% testing the function lmi_feas_tbx
% on the problem ||X+A||->min, l
% (A is rectangular, given, X is symmetric of size nxn)
%
% Written by ameg@mit.edu,  last modified October 13, 1997
if nargin<1, n=5; end   % default arguments
if nargin<2, cs=1; end

randn('state',cs);      % set the random generator seed
A=randn(n);             % generate the system

disp('*** Testing lmi software on an eigenvalue calculation problem')
disp('*')
disp(['*     answer should be = ' num2str(0.5*norm(A-A'))])
disp(' ')

abst_init_lmi;
x=symmetric(n);
[OO(n) A-x;A'-x OO(n)]<0;
lmi_mincx_tbx;

