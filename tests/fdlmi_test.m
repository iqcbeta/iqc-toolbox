% Test program for Frequency-Dependent LMIs solver
% 
% Fixed basis model approximation problem:
%
% Find a best H-inf approximation of G0=1/(s^3+2*s^2+s+1)
% using the basis 1, 1/(s+1), and 1/(s+2). 
% i.e.
% Let G=x0+x1/(s+1)+x2/(s+2). Find x0, x1, and x2
% such that || W*(G0 - G) ||_inf --> minimized, W is a weighting function
%
% written by cykao@mit.edu  on June 22 1999

%
% change call of "value" to "value_iqc" to avoid icompatibility with CONTROL MATLAB toolbox
%
% andres marcos, 27-October-2006 marcosa@aem.umn.edu



disp(' ')
disp(' approximate 1/(s^2+s+1) by 1, 1/(s+1), and 1/(s+0.7) ... ') 
disp(' ')

clear all
s=tf([1,0],1);
abst_init_fdlmi;

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

G0=1/(s*s+s+1);
G1=ss(-1,1,1,0);
G2=ss(-0.7,1,1,0);
x0=symmetric;
x1=symmetric;
x2=symmetric;
y =symmetric;
y>0;
Ga=x0+x1*G1+x2*G2;
H = (G0-Ga);
[y,H;H',y]>0;
fdlmi_mincx_tbx(y);

disp(' ')
disp(' the best approximation is: ')
value_iqc(Ga)
value_iqc(y)
norm(G0-value_iqc(Ga),inf)
bode(G0,value_iqc(Ga))


