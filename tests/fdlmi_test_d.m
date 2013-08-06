
disp(' ')
disp(' approximate 1/((z+0.6)*(z-0.3)) by 1, 1/(z+0.5), and 1/(z+0.7) ... ') 
disp(' ')

clear all
z=tf([1,0],1,-1);
abst_init_fdlmi(1);

% setlmioptions('lmilab')
% setlmioptions('yalmip','solver','sedumi')
% setlmioptions('yalmip','solver','sdpt3')

G0=1/((z+0.6)*(z-0.3));
G1=ss(-0.5,1,1,0,-1);
G2=ss(-0.7,1,1,0,-1);
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


