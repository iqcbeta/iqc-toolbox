% testing the "iqc_get_mlmi" functionality

clear all
s=tf([1 0],1);
G=(s-1)/((s+2)*(s+2));
abst_init_iqc;
lmitbx_options([0 0 0 0 1]);
w=signal;
v=G*w;
gain=iqc_gain_tbx(w,v);
[P,A,B,Sigma,MainLMI]=iqc_get_mlmi;

%% verifiying ... %%

[AG,BG,CG,DG]=ssdata(G);
T=[CG,DG;zeros(size(DG,2),size(CG,2)),eye(size(DG,2))];
Sigma_verify=T'*[1,0;0,-gain*gain]*T;

disp(' ')
disp('... the following should be very small number: ')
max(svd(Sigma_verify-Sigma))
