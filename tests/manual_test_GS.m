A0 = [0 0 1 0; 0 0 0 1; -1.25 1.25 0 0; 1.25 -1.25 0 0];
A1 = [0; 0; -sqrt(2)/2; sqrt(2)/2];
A2 = [0.75*sqrt(2) -0.75*sqrt(2) 0 0];
B1 = [0; 0; 1; 0];
B2 = [0; 0; 0; -1];
C  = [0 1 0 0];
Ac = [0 -0.7195 1 0; 0 -2.9732 0 1; 
      -2.5133 4.8548 -1.7287 -0.9616; 1.0063 -5.4097 -0.0081 0.0304];
Bc = [0.720; 2.973; -3.37; 4.419];
Cc = [-1.505 0.494 -1.738 -0.932];

Gp = ss(A0, eye(4), eye(4), zeros(4));
Gc = ss(Ac, Bc, Cc, 0);

abst_init_iqc
a = 0.1;
w1 = signal;
w2 = signal; 
w3 = signal;
u  = signal;
x  = Gp*(A1*w3+B1*(u+w1)+B2*w2);
v  = A2*x;
x2 = C*x;
uc = Gc*x2;
M  = variable;
M  > 0;
v'*(a^2*M)*v - w3'*M*w3 > 0;
u  == uc;
g  = iqc_gain_tbx([w1;w2],x2)