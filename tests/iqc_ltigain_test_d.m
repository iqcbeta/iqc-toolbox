
A0=[0 -0.5;1 1];
B0=[-6 0;1 0];
L0=[1 0];
delta=0.48;
st=1;

nA=size(A0,1);

abst_init_iqc(1)

setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sdpt3')

w=signal(size(B0,2));
G0=ss(A0,eye(nA),eye(nA),zeros(nA,nA),st);

p=signal;

x=G0*(B0*w+[0;delta]*p);
q=[0 1]*x;
v=L0*x;

p==iqc_ltigain(q);

iqc_gain_tbx(w,v)
