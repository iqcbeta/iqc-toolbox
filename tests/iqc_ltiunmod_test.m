% testing "iqc_ltiunmod": iqc description of an
% lti uncertain system with bounded norm.
%
% The L2 gain from f to v in the interconnection 
%
% v = M(f+w)   w = diag(delta1,delta2) v
%
% is compared for the cases in which each delta is
% (1)  a time variant uncertainty (using iqc_tvnorm)
% (2)  a linear time invariant uncertainty (using iqc_ltiunmod)
%
% Written by F.D'Amato,   last modified June 11, 1998

clear all
s=tf([1 0],1);
M=1/1.5*[s/(s+100)    0;
         s/(s+80)     0;
        10/(s+10)   8/(s+8)];

% case (1) --------------
abst_init_iqc

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

f=signal(2);
w1=signal;
w2=signal;

v=M*(f+[w1;w2]);
[waux1,x]=iqc_ltvnorm(v(1),1);
[waux2,x]=iqc_ltvnorm(v(2:3),1);
w1==waux1;
w2==waux2;

g=iqc_gain_tbx(f,v);
fprintf(1,' gain from f->v using iqc_ltvnorm : %f\n',g);

% case (2) --------------
abst_init_iqc
f=signal(2);
w1=signal;
w2=signal;
a=4;

v=M*(f+[w1;w2]);
[waux1,x]=iqc_ltiunmod(v(1),a,1);
[waux2,x]=iqc_ltiunmod(v(2:3),a,1);
w1==waux1;
w2==waux2;

g=iqc_gain_tbx(f,v);
fprintf(1,' gain from f->v using iqc_ltiunmod : %f\n',g);
