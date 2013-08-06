abst_init_iqc
% setlmioptions('yalmip','solver','sdpt3')
% setlmioptions('lmilab')

f=signal;
w=signal;
y=tf([1 1 7],[1 3.1 3.1 9])*w+f;
w==iqc_ltigain(y,[1,2.8]);
iqc_gain_tbx(f,y)
iqc_bode
