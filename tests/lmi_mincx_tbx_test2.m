function lmi_mincx_tbx_test2(n,cs)
% function lmi_mincx_tbx_test2(n,cs)
%
% testing the function lmi_mincx_tbx
% solving the problem x>a,x>b,trace(x)->min
% where a=a',b=b' are randomly generated n by n
%
% Written by ameg@mit.edu,  last modified October 13, 1997


%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu




if nargin<1, n=20; end   % default arguments
if nargin<2, cs=1; end

randn('state',cs);      % set the random generator seed
A=randn(n);             % generate the data
B=randn(n);

d=eig(0.5*(A'+A-B'-B));
answ=sum(d(find(d>0)))+trace(B);

disp('*** Testing lmi software on problem')
disp('*** x>a, x>b, trace(x)->min ')
disp('*')
disp(['*     answer = ' num2str(answ)])
disp(' ')

abst_init_lmi;

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

x=symmetric(n);
x>A;
x>B;
y=trace(x);
lmi_mincx_tbx(y);
X=value_iqc(x);
trace(X)
