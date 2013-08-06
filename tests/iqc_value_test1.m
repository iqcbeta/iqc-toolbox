% testing the "value" functionality

%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu



disp('The answer must be small ...')
s=tf([1 0],1);
G=(s-1)/((s+2)*(s+2));
H=(s-2)/((s+1)*(s+1));
abst_init_iqc;
lmitbx_options([0 0 0 0 1]);

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

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
norm(d,Inf)

